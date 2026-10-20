"""Download PubChem Compound SDF chunks from the NCBI HTTPS mirror.

PubChem publishes the full compound collection as ~350 gzipped SDF files of
up to 500,000 CIDs each under ``Compound/CURRENT-Full/SDF/``, each with an
``.md5`` sidecar. Downloads resume with HTTP ``Range`` requests and are
verified against the published MD5 before use.
"""

import hashlib
import logging
import os
import re
from collections.abc import Callable
from pathlib import Path

import httpx

from .. import __version__
from ..core.exceptions import DownloadError

logger = logging.getLogger(__name__)

DEFAULT_SDF_URL = "https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/SDF/"
CHUNK_RE = re.compile(r'href="(Compound_\d+_\d+\.sdf\.gz)"')
MD5_RE = re.compile(r"\b([0-9a-fA-F]{32})\b")
BLOCK = 1 << 20

ProgressFn = Callable[[str, int, int | None], None]
"""Called with (file name, bytes so far, total bytes or None)."""


class SdfMirror:
    """
    Client for the PubChem Compound SDF directory.

    Args:
        base_url: Directory URL ending in ``/``.
        client: Optional ``httpx.Client`` (tests, proxies).
        timeout: Read timeout in seconds.
    """

    def __init__(
        self,
        base_url: str = DEFAULT_SDF_URL,
        client: httpx.Client | None = None,
        timeout: float = 120.0,
    ) -> None:
        self.base_url = base_url if base_url.endswith("/") else base_url + "/"
        self._owns_client = client is None
        self._client = client or httpx.Client(
            timeout=httpx.Timeout(timeout, connect=30.0),
            headers={"User-Agent": f"pubchem-finder/{__version__}"},
            follow_redirects=True,
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
        try:
            response = self._client.get(self.base_url)
            response.raise_for_status()
        except httpx.HTTPError as exc:
            raise DownloadError(f"Could not list {self.base_url}: {exc}") from exc
        names = sorted(set(CHUNK_RE.findall(response.text)))
        if not names:
            raise DownloadError(f"No Compound_*.sdf.gz files found at {self.base_url}")
        return names

    def expected_md5(self, name: str) -> str:
        """The published MD5 of ``name``."""
        try:
            response = self._client.get(self.base_url + name + ".md5")
            response.raise_for_status()
        except httpx.HTTPError as exc:
            raise DownloadError(f"Could not fetch checksum for {name}: {exc}") from exc
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

        An existing verified file is reused; a partial ``.part`` file is resumed.

        Raises:
            DownloadError: network failure or checksum mismatch.
        """
        dest_dir.mkdir(parents=True, exist_ok=True)
        target = dest_dir / name
        md5 = self.expected_md5(name)
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
                self._fetch(name, part, progress)
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

    def _fetch(self, name: str, part: Path, progress: ProgressFn | None) -> None:
        offset = part.stat().st_size if part.exists() else 0
        headers = {"Range": f"bytes={offset}-"} if offset else {}
        url = self.base_url + name
        try:
            with self._client.stream("GET", url, headers=headers) as response:
                if response.status_code == 416:
                    raise _RangeNotSatisfiable
                if response.status_code not in (200, 206):
                    raise DownloadError(
                        f"GET {url} returned HTTP {response.status_code}"
                    )
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
            raise DownloadError(f"Download of {name} failed: {exc}") from exc


class _RangeNotSatisfiable(Exception):
    pass


def file_md5(path: Path) -> str:
    digest = hashlib.md5(usedforsecurity=False)
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(BLOCK), b""):
            digest.update(block)
    return digest.hexdigest()
