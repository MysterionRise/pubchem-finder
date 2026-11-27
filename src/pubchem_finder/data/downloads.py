"""Download PubChem Compound SDF chunks from the NCBI HTTPS mirror.

PubChem publishes the full compound collection as ~350 gzipped SDF files of
up to 500,000 CIDs each under ``Compound/CURRENT-Full/SDF/``, each with an
``.md5`` sidecar. Downloads resume with HTTP ``Range`` requests and are
verified against the published MD5 before use. Transient failures (a dropped
connection, a 503) are retried after a pause instead of ending a build that
may already have run for hours; the retry resumes from the bytes on disk.
"""

import hashlib
import logging
import os
import re
import time
from collections.abc import Callable
from pathlib import Path
from typing import TypeVar

import httpx

from .. import __version__
from ..core.exceptions import DownloadError, InvalidInput
from .pug_rest import no_cookies

logger = logging.getLogger(__name__)

DEFAULT_SDF_URL = "https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/SDF/"
CHUNK_RE = re.compile(r'href="(Compound_\d+_\d+\.sdf\.gz)"')
MD5_RE = re.compile(r"\b([0-9a-fA-F]{32})\b")
BLOCK = 1 << 20

DEFAULT_ATTEMPTS = 4
"""Attempts per file, the first one included."""
DEFAULT_RETRY_DELAY = 5.0
"""Seconds before the second attempt; every further pause is twice as long."""
TRANSIENT_STATUS = frozenset({408, 425, 429, 500, 502, 503, 504})
"""HTTP statuses a later attempt may not hit."""

T = TypeVar("T")

ProgressFn = Callable[[str, int, int | None], None]
"""Called with (file name, bytes so far, total bytes or None)."""


class SdfMirror:
    """
    Client for the PubChem Compound SDF directory.

    Args:
        base_url: Directory URL ending in ``/``.
        client: Optional ``httpx.Client`` (tests, proxies).
        timeout: Read timeout in seconds.
        attempts: Tries per file before a transient failure is reported, at least 1.
        retry_delay: Pause before the second try, doubled for each later one.
        sleep: The pause, replaced in tests to keep the suite instant.

    Raises:
        InvalidInput: ``attempts`` is below 1.
    """

    def __init__(
        self,
        base_url: str = DEFAULT_SDF_URL,
        client: httpx.Client | None = None,
        timeout: float = 120.0,
        *,
        attempts: int = DEFAULT_ATTEMPTS,
        retry_delay: float = DEFAULT_RETRY_DELAY,
        sleep: Callable[[float], None] = time.sleep,
    ) -> None:
        if attempts < 1:
            raise InvalidInput(f"attempts must be at least 1, not {attempts}")
        self.base_url = base_url if base_url.endswith("/") else base_url + "/"
        self.attempts = attempts
        self.retry_delay = retry_delay
        self._sleep = sleep
        self._owns_client = client is None
        self._client = client or httpx.Client(
            timeout=httpx.Timeout(timeout, connect=30.0),
            headers={"User-Agent": f"pubchem-finder/{__version__}"},
            follow_redirects=True,
            cookies=no_cookies(),
        )

    def close(self) -> None:
        if self._owns_client:
            self._client.close()

    def __enter__(self) -> "SdfMirror":
        return self

    def __exit__(self, *exc: object) -> None:
        self.close()

    def list_chunks(self) -> list[str]:
        """Names of all SDF chunks in the directory, in CID order."""
        response = self._retry(self._listing_page)
        names = sorted(set(CHUNK_RE.findall(response.text)))
        if not names:
            raise DownloadError(f"No Compound_*.sdf.gz files found at {self.base_url}")
        return names

    def _listing_page(self) -> httpx.Response:
        """The directory page, a transient failure being retryable like the rest."""
        try:
            response = self._client.get(self.base_url)
            response.raise_for_status()
        except httpx.HTTPError as exc:
            raise _failure(f"Could not list {self.base_url}: {exc}", exc) from exc
        return response

    def expected_md5(self, name: str) -> str:
        """The published MD5 of ``name``."""
        try:
            response = self._client.get(self.base_url + name + ".md5")
            response.raise_for_status()
        except httpx.HTTPError as exc:
            raise _failure(f"Could not fetch checksum for {name}: {exc}", exc) from exc
        match = MD5_RE.search(response.text)
        if not match:
            raise DownloadError(
                f"Unrecognized checksum file for {name}: {response.text!r}"
            )
        return match.group(1).lower()

    def download(
        self, name: str, dest_dir: Path, progress: ProgressFn | None = None
    ) -> Path:
        """
        Download ``name`` into ``dest_dir`` and verify its MD5.

        An existing verified file is reused; a partial ``.part`` file is resumed,
        also across retries: every attempt continues at the bytes already on disk.

        Raises:
            DownloadError: network failure or checksum mismatch.
        """
        dest_dir.mkdir(parents=True, exist_ok=True)
        target = dest_dir / name
        md5 = self._retry(lambda: self.expected_md5(name))
        if target.exists():
            if file_md5(target) == md5:
                return target
            logger.warning(
                "%s exists but its checksum is wrong; downloading again", name
            )
            target.unlink()

        part = dest_dir / (name + ".part")
        for attempt in (1, 2):
            try:
                self._retry(lambda: self._fetch(name, part, progress))
                break
            except _RangeNotSatisfiable:
                part.unlink(missing_ok=True)
                if attempt == 2:
                    raise DownloadError(f"Server refused to resume {name}") from None
        actual = file_md5(part)
        if actual != md5:
            part.unlink(missing_ok=True)
            raise DownloadError(
                f"Checksum mismatch for {name}: expected {md5}, got {actual}"
            )
        os.replace(part, target)
        return target

    def _retry(self, fetch: Callable[[], T]) -> T:
        """
        Run ``fetch``, pausing and trying again after a transient failure.

        A long build must not die because one request timed out, while a
        permanent answer (a 404, a checksum file we cannot read) is final on the
        first try: retrying it only spends the build's time.

        Raises:
            DownloadError: the last failure, once every attempt has been used.
        """
        delay = self.retry_delay
        for attempt in range(1, self.attempts + 1):
            try:
                return fetch()
            except _TransientError as exc:
                if attempt == self.attempts:
                    raise DownloadError(
                        f"{exc}; gave up after {attempt} attempts"
                    ) from exc
                logger.warning(
                    "%s; retrying in %gs (attempt %d of %d)",
                    exc,
                    delay,
                    attempt,
                    self.attempts,
                )
                self._sleep(delay)
                delay *= 2
        raise AssertionError("unreachable")  # pragma: no cover

    def _fetch(self, name: str, part: Path, progress: ProgressFn | None) -> None:
        offset = part.stat().st_size if part.exists() else 0
        headers = {"Range": f"bytes={offset}-"} if offset else {}
        url = self.base_url + name
        try:
            with self._client.stream("GET", url, headers=headers) as response:
                if response.status_code == 416:
                    raise _RangeNotSatisfiable
                if response.status_code not in (200, 206):
                    raise _status_failure(url, response.status_code)
                if response.status_code == 200:
                    offset = 0  # server ignored the Range header: start over
                length = response.headers.get("Content-Length")
                total = offset + int(length) if length else None
                done = offset
                with open(part, "ab" if offset else "wb") as fh:
                    for block in response.iter_bytes(BLOCK):
                        fh.write(block)
                        done += len(block)
                        if progress:
                            progress(name, done, total)
        except httpx.HTTPError as exc:
            raise _failure(f"Download of {name} failed: {exc}", exc) from exc


class _RangeNotSatisfiable(Exception):
    pass


class _TransientError(DownloadError):
    """A failure a later attempt may not hit: a dropped connection, a timeout, or a
    408, 425, 429, 500, 502, 503 or 504 answer."""


def _failure(message: str, exc: httpx.HTTPError) -> DownloadError:
    """Wrap ``exc``, keeping network blips and server hiccups retryable."""
    if isinstance(exc, httpx.TransportError):
        return _TransientError(message)
    response = getattr(exc, "response", None)
    if getattr(response, "status_code", None) in TRANSIENT_STATUS:
        return _TransientError(message)
    return DownloadError(message)


def _status_failure(url: str, status: int) -> DownloadError:
    """The error for a download answered with ``status``."""
    message = f"GET {url} returned HTTP {status}"
    return (
        _TransientError(message)
        if status in TRANSIENT_STATUS
        else DownloadError(message)
    )


def file_md5(path: Path) -> str:
    digest = hashlib.md5(usedforsecurity=False)
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(BLOCK), b""):
            digest.update(block)
    return digest.hexdigest()
