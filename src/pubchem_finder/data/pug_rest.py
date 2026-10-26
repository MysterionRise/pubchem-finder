"""Client for PubChem's PUG REST and PUG View web services.

PubChem asks for no more than 5 requests per second and returns HTTP 503
("ServerBusy") when it is overloaded; this client throttles and retries.
Successful responses can be cached on disk (see :mod:`.cache`).
"""

import json
import logging
import threading
import time
from collections.abc import Mapping, Sequence
from typing import Any, Literal, get_args
from urllib.parse import quote

import httpx

from .. import __version__
from ..core.exceptions import PubChemFinderError
from ..models import CompoundRecord
from .cache import ResponseCache

logger = logging.getLogger(__name__)

DEFAULT_BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
DEFAULT_VIEW_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view"

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

IdentityType = Literal[
    "same_stereo_isotope",
    "same_stereo",
    "same_isotope",
    "same_connectivity",
    "same_tautomer",
    "nonconflicting_stereo",
]
IDENTITY_TYPES: tuple[str, ...] = get_args(IdentityType)

XREF_TYPES: tuple[str, ...] = ("RegistryID", "PubMedID", "PatentID")

# Compound page sections offered to the assistant. These are PubChem TOC headings;
# requesting a parent heading returns all of its subsections.
AnnotationHeading = Literal[
    "GHS Classification",
    "Hazards Identification",
    "Pharmacology and Biochemistry",
    "Mechanism of Action",
    "Drug and Medication Information",
    "Experimental Properties",
    "Toxicity",
    "Use and Manufacturing",
    "Food Additives and Ingredients",
]
ANNOTATION_HEADINGS: tuple[str, ...] = get_args(AnnotationHeading)

RETRY_STATUSES = frozenset({429, 503})
MAX_SEARCH_RECORDS = 200


class PubChemError(PubChemFinderError):
    """PubChem returned an error or an unexpected response."""


class CompoundNotFound(PubChemError):
    """No compound matches the query."""


class BadQuery(PubChemError):
    """PubChem rejected the query (malformed identifier)."""


class PubChemUnavailable(PubChemError):
    """PubChem is busy or down and retries were exhausted."""


class PubChemTimeout(PubChemError):
    """A PubChem search ran past the server's time limit."""


class AnnotationNotFound(PubChemError):
    """The compound page has no section with the requested heading."""

    def __init__(self, message: str, available: Sequence[str] = ()) -> None:
        super().__init__(message)
        self.available = list(available)


class PubChemClient:
    """
    Synchronous PUG REST client.

    Args:
        base_url: PUG REST root.
        timeout: Per-request timeout in seconds.
        client: Optional preconfigured ``httpx.Client`` (tests, custom proxies).
        min_interval: Minimum seconds between requests (PubChem allows 5/s).
        max_retries: Attempts on 429/503 before raising ``PubChemUnavailable``.
        cache: Optional on-disk cache for successful responses.
        view_url: PUG View root.
    """

    def __init__(
        self,
        base_url: str = DEFAULT_BASE_URL,
        timeout: float = 30.0,
        client: httpx.Client | None = None,
        min_interval: float = 0.2,
        max_retries: int = 3,
        cache: ResponseCache | None = None,
        view_url: str = DEFAULT_VIEW_URL,
    ) -> None:
        self.base_url = base_url.rstrip("/")
        self.view_url = view_url.rstrip("/")
        self.cache = cache
        self.min_interval = min_interval
        self.max_retries = max_retries
        self._owns_client = client is None
        self._client = client or httpx.Client(
            timeout=timeout,
            headers={"User-Agent": f"pubchem-finder/{__version__}"},
        )
        self._last_request: float | None = None
        # MCP tools run on worker threads and share one client, so the throttle
        # must serialize them; the wait happens inside the lock on purpose.
        self._throttle_lock = threading.Lock()

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
            payload = self.request("POST", url, data={"smiles": query})
        else:
            url = f"{self.base_url}/compound/{by}/{quote(query, safe='')}/property/{props}/JSON"
            payload = self.request("GET", url)

        rows = _property_rows(payload)
        if not rows:
            raise CompoundNotFound(f"No PubChem compound found for {by}={query!r}")
        return [_to_record(row) for row in rows[:limit]]

    # ------------------------------------------------------------------ #
    # Structure searches over all of PubChem
    # ------------------------------------------------------------------ #

    def structure_search(
        self,
        kind: Literal["fastsimilarity_2d", "fastsubstructure", "fastidentity"],
        smiles: str,
        max_records: int = 25,
        **params: Any,
    ) -> list[CompoundRecord]:
        """
        Run one of PubChem's synchronous ("fast") structure searches and return
        property records for the hits, in PubChem's order. No hits -> ``[]``.

        Raises:
            BadQuery: empty SMILES or out-of-range ``max_records``.
            PubChemTimeout: the search hit PubChem's time limit.
        """
        smiles = smiles.strip()
        if not smiles:
            raise BadQuery("SMILES must not be empty")
        if not 1 <= max_records <= MAX_SEARCH_RECORDS:
            raise BadQuery(f"max_records must be between 1 and {MAX_SEARCH_RECORDS}")
        props = ",".join(PROPERTIES)
        url = f"{self.base_url}/compound/{kind}/smiles/property/{props}/JSON"
        query = {**params, "MaxRecords": max_records}
        try:
            payload = self.request("POST", url, params=query, data={"smiles": smiles})
        except CompoundNotFound:
            return []
        return [_to_record(row) for row in _property_rows(payload)[:max_records]]

    def similarity_search(
        self, smiles: str, threshold: int = 90, max_records: int = 25
    ) -> list[CompoundRecord]:
        """Compounds with PubChem fingerprint Tanimoto >= ``threshold`` percent."""
        if not 0 <= threshold <= 100:
            raise BadQuery("threshold must be between 0 and 100")
        return self.structure_search(
            "fastsimilarity_2d", smiles, max_records, Threshold=threshold
        )

    def substructure_search(
        self, smiles: str, max_records: int = 25
    ) -> list[CompoundRecord]:
        """Compounds that contain ``smiles`` as a substructure."""
        return self.structure_search("fastsubstructure", smiles, max_records)

    def exact_match(
        self,
        smiles: str,
        identity_type: IdentityType = "same_stereo_isotope",
        max_records: int = 25,
    ) -> list[CompoundRecord]:
        """Compounds identical to ``smiles`` under ``identity_type``."""
        if identity_type not in IDENTITY_TYPES:
            raise BadQuery(f"identity_type must be one of {IDENTITY_TYPES}")
        return self.structure_search(
            "fastidentity", smiles, max_records, identity_type=identity_type
        )

    # ------------------------------------------------------------------ #
    # Per-compound information
    # ------------------------------------------------------------------ #

    def compound_json(self, cid: int, operation: str) -> dict[str, Any]:
        """``GET /compound/cid/{cid}/{operation}/JSON``."""
        if cid < 1:
            raise BadQuery("CID must be a positive integer")
        return self.request(
            "GET", f"{self.base_url}/compound/cid/{cid}/{operation}/JSON"
        )

    def view_data(self, cid: int, heading: str) -> dict[str, Any]:
        """PUG View record for one heading of a compound page."""
        if cid < 1:
            raise BadQuery("CID must be a positive integer")
        return self.request(
            "GET",
            f"{self.view_url}/data/compound/{cid}/JSON",
            params={"heading": heading},
        )

    def view_index(self, cid: int) -> dict[str, Any]:
        """PUG View table of contents of a compound page."""
        if cid < 1:
            raise BadQuery("CID must be a positive integer")
        return self.request("GET", f"{self.view_url}/index/compound/{cid}/JSON")

    # ------------------------------------------------------------------ #
    # Transport
    # ------------------------------------------------------------------ #

    def request(
        self,
        method: str,
        url: str,
        params: Mapping[str, Any] | None = None,
        data: Mapping[str, Any] | None = None,
    ) -> dict[str, Any]:
        """Send a request (or answer it from the cache) and return the JSON body."""
        key = ResponseCache.key(method, url, params, data) if self.cache else ""
        if self.cache is not None:
            cached = self.cache.get(key)
            if cached is not None:
                result: dict[str, Any] = json.loads(cached)
                return result
        response = self._send(method, url, params, data)
        payload = self._handle(response)
        if self.cache is not None:
            self.cache.put(key, response.content)
        return payload

    def _send(
        self,
        method: str,
        url: str,
        params: Mapping[str, Any] | None,
        data: Mapping[str, Any] | None,
    ) -> httpx.Response:
        delay = 0.5
        for attempt in range(1, self.max_retries + 1):
            self._throttle()
            try:
                response = self._client.request(
                    method,
                    url,
                    params=dict(params) if params else None,
                    data=dict(data) if data else None,
                )
            except httpx.TimeoutException as exc:
                raise PubChemTimeout(f"PubChem did not answer in time: {exc}") from exc
            except httpx.HTTPError as exc:
                raise PubChemError(f"Request to PubChem failed: {exc}") from exc

            if (
                response.status_code in RETRY_STATUSES
                and attempt < self.max_retries
                and not _is_timeout(response)
            ):
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
            return response
        raise AssertionError("unreachable")  # pragma: no cover

    def _throttle(self) -> None:
        with self._throttle_lock:
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
            try:
                data: dict[str, Any] = response.json()
            except ValueError as exc:
                raise PubChemError(f"PubChem returned invalid JSON: {exc}") from exc
            return data
        message = _fault_message(response)
        if status == 504 or _is_timeout(response):
            raise PubChemTimeout(f"PubChem search timed out: {message}")
        if status == 404:
            raise CompoundNotFound(message)
        if status == 400:
            raise BadQuery(message)
        if status in RETRY_STATUSES:
            raise PubChemUnavailable(f"PubChem is busy ({status}): {message}")
        raise PubChemError(f"PubChem returned HTTP {status}: {message}")


def _property_rows(payload: dict[str, Any]) -> list[dict[str, Any]]:
    try:
        rows: list[dict[str, Any]] = list(payload["PropertyTable"]["Properties"])
    except (KeyError, TypeError) as exc:
        raise PubChemError(f"Unexpected PubChem response shape: {exc}") from exc
    return rows


def _is_timeout(response: httpx.Response) -> bool:
    try:
        code = str(response.json()["Fault"].get("Code", ""))
    except (ValueError, KeyError, TypeError, AttributeError):
        return False
    return "timeout" in code.lower()


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
