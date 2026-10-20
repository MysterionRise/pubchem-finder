"""Minimal client for PubChem's PUG REST web service.

Only the compound property table is used:
``/compound/{namespace}/{identifier}/property/{properties}/JSON``.

PubChem asks for no more than 5 requests per second and returns HTTP 503
("ServerBusy") when it is overloaded; this client throttles and retries.
"""

import logging
import time
from collections.abc import Sequence
from typing import Any, Literal, get_args
from urllib.parse import quote

import httpx

from .. import __version__
from ..core.exceptions import PubChemFinderError
from ..models import CompoundRecord

logger = logging.getLogger(__name__)

DEFAULT_BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

# PubChem renamed CanonicalSMILES -> ConnectivitySMILES and IsomericSMILES -> SMILES
# in 2025; the old names are still accepted as aliases but echoed under the new ones.
PROPERTIES: tuple[str, ...] = (
    "Title",
    "IUPACName",
    "SMILES",
    "ConnectivitySMILES",
    "InChIKey",
    "MolecularFormula",
    "MolecularWeight",
    "XLogP",
    "HBondDonorCount",
    "HBondAcceptorCount",
    "TPSA",
)

# PUG REST property name -> CompoundRecord field
FIELD_MAP: dict[str, str] = {
    "CID": "cid",
    "Title": "title",
    "IUPACName": "iupac_name",
    "SMILES": "smiles",
    "ConnectivitySMILES": "connectivity_smiles",
    "InChIKey": "inchikey",
    "MolecularFormula": "molecular_formula",
    "MolecularWeight": "molecular_weight",
    "XLogP": "xlogp",
    "HBondDonorCount": "h_bond_donors",
    "HBondAcceptorCount": "h_bond_acceptors",
    "TPSA": "tpsa",
}

Namespace = Literal["name", "cid", "smiles", "inchikey"]
NAMESPACES: tuple[str, ...] = get_args(Namespace)

RETRY_STATUSES = frozenset({429, 503})


class PubChemError(PubChemFinderError):
    """PubChem returned an error or an unexpected response."""


class CompoundNotFound(PubChemError):
    """No compound matches the query."""


class BadQuery(PubChemError):
    """PubChem rejected the query (malformed identifier)."""


class PubChemUnavailable(PubChemError):
    """PubChem is busy or down and retries were exhausted."""


class PubChemClient:
    """
    Synchronous PUG REST client.

    Args:
        base_url: PUG REST root.
        timeout: Per-request timeout in seconds.
        client: Optional preconfigured ``httpx.Client`` (tests, custom proxies).
        min_interval: Minimum seconds between requests (PubChem allows 5/s).
        max_retries: Attempts on 429/503 before raising ``PubChemUnavailable``.
    """

    def __init__(
        self,
        base_url: str = DEFAULT_BASE_URL,
        timeout: float = 20.0,
        client: httpx.Client | None = None,
        min_interval: float = 0.2,
        max_retries: int = 3,
    ) -> None:
        self.base_url = base_url.rstrip("/")
        self.min_interval = min_interval
        self.max_retries = max_retries
        self._owns_client = client is None
        self._client = client or httpx.Client(
            timeout=timeout,
            headers={"User-Agent": f"pubchem-finder/{__version__}"},
        )
        self._last_request: float | None = None

    def close(self) -> None:
        if self._owns_client:
            self._client.close()

    def __enter__(self) -> "PubChemClient":
        return self

    def __exit__(self, *exc: object) -> None:
        self.close()

    def lookup(
        self, query: str, by: Namespace = "name", limit: int = 5
    ) -> list[CompoundRecord]:
        """
        Look up compounds and return their property records.

        Args:
            query: Name, CID (or comma-separated CIDs), SMILES or InChIKey.
            by: Which identifier namespace ``query`` is in.
            limit: Maximum number of records to return.

        Raises:
            BadQuery: empty query or unknown namespace, or PubChem rejected it.
            CompoundNotFound: nothing matched.
            PubChemUnavailable: PubChem kept answering 429/503.
            PubChemError: any other failure.
        """
        query = query.strip()
        if not query:
            raise BadQuery("Query must not be empty")
        if by not in NAMESPACES:
            raise BadQuery(f"Unknown namespace {by!r}; use one of {NAMESPACES}")
        if limit < 1:
            raise BadQuery("limit must be >= 1")

        props = ",".join(PROPERTIES)
        if by == "smiles":
            # SMILES can contain '/' and '#', which do not survive in a URL path.
            url = f"{self.base_url}/compound/smiles/property/{props}/JSON"
            payload = self._request("POST", url, data={"smiles": query})
        else:
            url = f"{self.base_url}/compound/{by}/{quote(query, safe='')}/property/{props}/JSON"
            payload = self._request("GET", url)

        try:
            rows: Sequence[dict[str, Any]] = payload["PropertyTable"]["Properties"]
        except (KeyError, TypeError) as exc:
            raise PubChemError(f"Unexpected PubChem response shape: {exc}") from exc
        if not rows:
            raise CompoundNotFound(f"No PubChem compound found for {by}={query!r}")
        return [_to_record(row) for row in rows[:limit]]

    def _request(self, method: str, url: str, **kwargs: Any) -> dict[str, Any]:
        delay = 0.5
        for attempt in range(1, self.max_retries + 1):
            self._throttle()
            try:
                response = self._client.request(method, url, **kwargs)
            except httpx.HTTPError as exc:
                raise PubChemError(f"Request to PubChem failed: {exc}") from exc

            if response.status_code in RETRY_STATUSES and attempt < self.max_retries:
                logger.warning(
                    "PubChem answered %s (attempt %d/%d), retrying in %.1fs",
                    response.status_code,
                    attempt,
                    self.max_retries,
                    delay,
                )
                time.sleep(delay)
                delay *= 2
                continue
            return self._handle(response)
        raise AssertionError("unreachable")  # pragma: no cover

    def _throttle(self) -> None:
        now = time.monotonic()
        if self._last_request is not None:
            wait = self.min_interval - (now - self._last_request)
            if wait > 0:
                time.sleep(wait)
        self._last_request = time.monotonic()

    @staticmethod
    def _handle(response: httpx.Response) -> dict[str, Any]:
        status = response.status_code
        if status == 200:
            data: dict[str, Any] = response.json()
            return data
        message = _fault_message(response)
        if status == 404:
            raise CompoundNotFound(message)
        if status == 400:
            raise BadQuery(message)
        if status in RETRY_STATUSES:
            raise PubChemUnavailable(f"PubChem is busy ({status}): {message}")
        raise PubChemError(f"PubChem returned HTTP {status}: {message}")


def _to_record(row: dict[str, Any]) -> CompoundRecord:
    data = {FIELD_MAP[key]: value for key, value in row.items() if key in FIELD_MAP}
    try:
        return CompoundRecord.model_validate(data)
    except ValueError as exc:
        raise PubChemError(f"Unexpected PubChem property values: {exc}") from exc


def _fault_message(response: httpx.Response) -> str:
    try:
        fault = response.json()["Fault"]
    except (ValueError, KeyError, TypeError):
        return response.text.strip()[:200] or f"HTTP {response.status_code}"
    parts = [str(fault.get("Message", "")).strip()]
    parts.extend(str(d) for d in fault.get("Details", []) if d)
    return "; ".join(p for p in parts if p) or str(fault.get("Code", "unknown fault"))
