"""Framework-free operations shared by the CLI and the MCP server.

Every function here takes plain values, returns a Pydantic model and raises
``PubChemFinderError`` subclasses. The CLI and the MCP server only translate
those into their own output and error formats.
"""

import logging
import sqlite3
import threading
from pathlib import Path
from typing import Literal

import rdkit
from rdkit import Chem, DataStructs, rdBase

from . import drawing, medchem
from .chemistry import MoleculeStandardizer, calculate_properties, is_druglike
from .config import cache_enabled, cache_path, resolve_home
from .core.exceptions import MoleculeValidationError, PubChemFinderError
from .data import parsers
from .data.cache import CacheStats, ResponseCache
from .data.pug_rest import (
    XREF_TYPES,
    AnnotationHeading,
    AnnotationNotFound,
    BadQuery,
    CompoundNotFound,
    IdentityType,
    Namespace,
    PubChemClient,
    PubChemError,
)
from .fingerprints import MACCSFingerprint, MorganFingerprint
from .index.engine import LocalIndex, Metric
from .models import (
    AnnotationSection,
    BioactivitySummary,
    Comparison,
    CompoundRecord,
    CrossRefs,
    Description,
    IndexedCompound,
    IndexStatus,
    LookupResult,
    MedchemProfile,
    Properties,
    RemoteHit,
    RemoteSearchResults,
    SimilarityHit,
    SimilarityResults,
    Standardized,
    SubstructureResults,
    Synonyms,
)

logger = logging.getLogger(__name__)

# RDKit writes parse warnings to stderr; keep server output quiet.
rdBase.DisableLog("rdApp.*")

_ECFP4 = MorganFingerprint(radius=2, n_bits=2048)
_MACCS = MACCSFingerprint()
_STANDARDIZER = MoleculeStandardizer()

QueryType = Literal["smiles", "smarts"]

_client_lock = threading.Lock()
_client: PubChemClient | None = None
_indexes: dict[Path, LocalIndex] = {}


def get_client() -> PubChemClient:
    """Process-wide PubChem client, created on first use, with the on-disk
    response cache unless ``PUBCHEM_FINDER_CACHE=0``."""
    global _client
    with _client_lock:
        if _client is None:
            cache = None
            if cache_enabled():
                try:
                    cache = ResponseCache(cache_path())
                except (OSError, sqlite3.Error) as exc:
                    logger.warning("PubChem response cache disabled: %s", exc)
            _client = PubChemClient(cache=cache)
        return _client


def parse_smiles(smiles: str) -> Chem.Mol:
    """Parse SMILES or raise ``MoleculeValidationError``."""
    mol: Chem.Mol | None = (
        Chem.MolFromSmiles(smiles.strip()) if smiles.strip() else None
    )
    if mol is None:
        raise MoleculeValidationError(f"Invalid SMILES: {smiles!r}")
    return mol


def parse_query(query: str, query_type: QueryType = "smiles") -> Chem.Mol:
    """Parse a substructure query given as SMILES or SMARTS."""
    if query_type == "smiles":
        return parse_smiles(query)
    if query_type != "smarts":
        raise PubChemFinderError(f"Unknown query_type {query_type!r}")
    mol: Chem.Mol | None = Chem.MolFromSmarts(query.strip()) if query.strip() else None
    if mol is None:
        raise MoleculeValidationError(f"Invalid SMARTS: {query!r}")
    # Pattern fingerprints need ring information, which SMARTS parsing skips.
    mol.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(mol)
    return mol


def lookup(query: str, by: Namespace = "name", limit: int = 5) -> LookupResult:
    """Look compounds up in PubChem."""
    compounds = get_client().lookup(query, by=by, limit=limit)
    return LookupResult(query=query, by=by, compounds=compounds)


def properties(smiles: str) -> Properties:
    """RDKit descriptors and Lipinski check for one structure."""
    mol = parse_smiles(smiles)
    props = calculate_properties(mol)
    return Properties.model_validate(
        {
            "input_smiles": smiles,
            "canonical_smiles": Chem.MolToSmiles(mol),
            "lipinski_ok": is_druglike(mol),
            **props,
        }
    )


def standardize(smiles: str) -> Standardized:
    """Standardize a structure (cleanup, parent fragment, neutralize)."""
    mol = parse_smiles(smiles)
    std = _STANDARDIZER.standardize(mol)
    if std is None:
        raise MoleculeValidationError(f"Could not standardize {smiles!r}")
    result = str(Chem.MolToSmiles(std))
    return Standardized(
        input_smiles=smiles,
        smiles=result,
        changed=result != Chem.MolToSmiles(mol),
    )


def compare(smiles_a: str, smiles_b: str) -> Comparison:
    """Tanimoto similarity between two structures on ECFP4 and MACCS."""
    mol_a = parse_smiles(smiles_a)
    mol_b = parse_smiles(smiles_b)
    return Comparison(
        smiles_a=smiles_a,
        smiles_b=smiles_b,
        ecfp4_tanimoto=float(
            DataStructs.TanimotoSimilarity(
                _ECFP4.generate_bitvect(mol_a), _ECFP4.generate_bitvect(mol_b)
            )
        ),
        maccs_tanimoto=float(
            DataStructs.TanimotoSimilarity(
                _MACCS.generate_bitvect(mol_a), _MACCS.generate_bitvect(mol_b)
            )
        ),
    )


def get_index(home: Path | str | None = None) -> LocalIndex:
    """The local index for ``home`` (default: configured data directory), cached."""
    path = resolve_home(home)
    with _client_lock:
        if path not in _indexes:
            _indexes[path] = LocalIndex(path)
        return _indexes[path]


def similarity_search(
    smiles: str,
    threshold: float = 0.7,
    top_k: int = 20,
    metric: Metric = "tanimoto",
    home: Path | str | None = None,
    group_by_scaffold: bool = False,
) -> SimilarityResults:
    """Most similar compounds in the local index."""
    if not 0.0 < threshold <= 1.0:
        raise PubChemFinderError("threshold must be in (0, 1]")
    if not 1 <= top_k <= 1000:
        raise PubChemFinderError("top_k must be between 1 and 1000")
    if metric not in ("tanimoto", "dice", "cosine"):
        raise PubChemFinderError(f"Unknown metric {metric!r}")
    mol = parse_smiles(smiles)
    index = get_index(home)
    manifest = index.manifest
    pairs = index.search(mol, threshold=threshold, top_k=top_k, metric=metric)
    meta = index.metadata([cid for cid, _ in pairs])
    hits = [
        SimilarityHit.model_validate(
            {**meta.get(cid, {}), "cid": cid, "similarity": sim}
        )
        for cid, sim in pairs
    ]
    return SimilarityResults(
        query_smiles=smiles,
        threshold=threshold,
        metric=metric,
        fingerprint=manifest.fingerprint,
        n_indexed=manifest.n_compounds,
        hits=hits,
        scaffolds=medchem.group_by_scaffold(hits) if group_by_scaffold else None,
    )


def substructure_search(
    query: str,
    query_type: QueryType = "smiles",
    max_hits: int = 50,
    max_candidates: int = 100_000,
    home: Path | str | None = None,
) -> SubstructureResults:
    """Compounds in the local index that contain ``query``, lowest CID first."""
    if not 1 <= max_hits <= 1000:
        raise PubChemFinderError("max_hits must be between 1 and 1000")
    if max_candidates < 1:
        raise PubChemFinderError("max_candidates must be at least 1")
    mol = parse_query(query, query_type)
    index = get_index(home)
    found = index.substructure(mol, max_hits=max_hits, max_candidates=max_candidates)
    sub = index.manifest.substructure
    warning = None
    if sub is not None and sub.rdkit_version != rdkit.__version__:
        warning = (
            f"The substructure index was built with RDKit {sub.rdkit_version}, but "
            f"RDKit {rdkit.__version__} is running. Pattern fingerprints can differ "
            "between releases, so matches may be missing. Rebuild it with "
            "`pubchem-finder index build --with-substructure`."
        )
    return SubstructureResults(
        query=query,
        query_type=query_type,
        n_indexed=sub.n_compounds if sub is not None else 0,
        n_candidates=found.n_candidates,
        n_checked=found.n_checked,
        exhaustive=found.exhaustive,
        hits=[IndexedCompound.model_validate(row) for row in found.rows],
        warning=warning,
    )


def index_status(home: Path | str | None = None) -> IndexStatus:
    """Describe the local index without loading it."""
    index = get_index(home)
    paths = index.paths
    built = paths.completed_chunks()
    manifest = paths.load_manifest()
    if manifest is None:
        return IndexStatus(
            state="not_built",
            home=str(paths.home),
            chunks_built=len(built),
            hint=(
                "Build with `pubchem-finder index build --limit 1` for a ~500k-compound "
                "sample, or `--all` for all of PubChem."
            ),
        )
    sub = manifest.substructure if paths.pattern_fingerprints.exists() else None
    return IndexStatus(
        state="ready",
        home=str(paths.home),
        n_compounds=manifest.n_compounds,
        n_failed=manifest.n_failed,
        fingerprint=manifest.fingerprint,
        built_at=manifest.built_at.isoformat(timespec="seconds"),
        size_bytes=paths.fingerprints.stat().st_size,
        memory_bytes=manifest.memory_bytes(),
        chunks_in_index=len(manifest.chunks),
        chunks_built=len(built),
        loaded=index.loaded,
        mode=None
        if index.in_memory is None
        else ("memory" if index.in_memory else "disk"),
        substructure=sub is not None,
        substructure_size_bytes=paths.pattern_fingerprints.stat().st_size
        if sub is not None
        else 0,
        substructure_memory_bytes=sub.memory_bytes() if sub is not None else 0,
        hint=None
        if sub is not None
        else (
            "Substructure search is off; add it with "
            "`pubchem-finder index build --with-substructure` (no download)."
        ),
    )


# --------------------------------------------------------------------------- #
# PubChem knowledge tools
# --------------------------------------------------------------------------- #

SIMILARITY_NOTE = (
    "PubChem selected these compounds with its own 881-bit substructure fingerprint "
    "(Tanimoto >= {threshold}%), which does not report a score per hit. ecfp4_similarity "
    "is computed here on ECFP4, the same fingerprint as the local similarity_search; "
    "hits are sorted by it."
)
SUBSTRUCTURE_NOTE = (
    "Every hit contains the query as a substructure. PubChem returns at most max_records "
    "hits; they are sorted by ecfp4_similarity to the query, so close analogues come first."
)
IDENTITY_NOTE = "Compounds PubChem considers identical to the query under identity_type={identity_type}."


def _rank(query: Chem.Mol, records: list[CompoundRecord]) -> list[RemoteHit]:
    query_fp = _ECFP4.generate_bitvect(query)
    hits: list[RemoteHit] = []
    for record in records:
        similarity = None
        smiles = record.smiles or record.connectivity_smiles
        mol = Chem.MolFromSmiles(smiles) if smiles else None
        if mol is not None:
            similarity = round(
                float(
                    DataStructs.TanimotoSimilarity(
                        query_fp, _ECFP4.generate_bitvect(mol)
                    )
                ),
                4,
            )
        hits.append(
            RemoteHit(
                **record.model_dump(exclude={"pubchem_url"}),
                ecfp4_similarity=similarity,
            )
        )
    hits.sort(
        key=lambda h: -1.0 if h.ecfp4_similarity is None else h.ecfp4_similarity,
        reverse=True,
    )
    return hits


def pubchem_similarity_search(
    smiles: str, threshold: int = 90, max_records: int = 25
) -> RemoteSearchResults:
    """Similar compounds across all of PubChem (PubChem's fast 2D similarity search)."""
    mol = parse_smiles(smiles)
    records = get_client().similarity_search(
        smiles, threshold=threshold, max_records=max_records
    )
    return RemoteSearchResults(
        query_smiles=smiles,
        search="similarity",
        threshold=threshold,
        max_records=max_records,
        truncated=len(records) >= max_records,
        hits=_rank(mol, records),
        note=SIMILARITY_NOTE.format(threshold=threshold),
    )


def pubchem_substructure_search(
    smiles: str, max_records: int = 25
) -> RemoteSearchResults:
    """Compounds across all of PubChem that contain ``smiles``."""
    mol = parse_smiles(smiles)
    records = get_client().substructure_search(smiles, max_records=max_records)
    return RemoteSearchResults(
        query_smiles=smiles,
        search="substructure",
        max_records=max_records,
        truncated=len(records) >= max_records,
        hits=_rank(mol, records),
        note=SUBSTRUCTURE_NOTE,
    )


def pubchem_exact_match(
    smiles: str, identity_type: IdentityType = "same_stereo_isotope"
) -> RemoteSearchResults:
    """PubChem compounds identical to ``smiles``."""
    mol = parse_smiles(smiles)
    max_records = 25
    records = get_client().exact_match(
        smiles, identity_type=identity_type, max_records=max_records
    )
    return RemoteSearchResults(
        query_smiles=smiles,
        search="identity",
        max_records=max_records,
        truncated=len(records) >= max_records,
        hits=_rank(mol, records),
        note=IDENTITY_NOTE.format(identity_type=identity_type),
    )


def pubchem_synonyms(cid: int, limit: int = 25) -> Synonyms:
    """Names, CAS numbers and other identifiers PubChem lists for a CID."""
    if not 1 <= limit <= 500:
        raise BadQuery("limit must be between 1 and 500")
    payload = get_client().compound_json(cid, "synonyms")
    try:
        names = parsers.parse_synonyms(payload)
    except parsers.PayloadError as exc:
        raise PubChemError(f"Unexpected synonyms response: {exc}") from exc
    return Synonyms(
        cid=cid,
        synonyms=names[:limit],
        cas_numbers=[n for n in names if parsers.is_cas_number(n)][:10],
        total=len(names),
    )


def pubchem_xrefs(cid: int, sample: int = 10) -> CrossRefs:
    """Registry IDs, PubMed articles and patents linked to a CID."""
    if not 1 <= sample <= 100:
        raise BadQuery("sample must be between 1 and 100")
    try:
        payload = get_client().compound_json(cid, "xrefs/" + ",".join(XREF_TYPES))
        refs = parsers.parse_xrefs(payload, XREF_TYPES)
    except CompoundNotFound:
        refs = {t: [] for t in XREF_TYPES}
    except parsers.PayloadError as exc:
        raise PubChemError(f"Unexpected xrefs response: {exc}") from exc
    return CrossRefs(
        cid=cid,
        counts={t: len(v) for t, v in refs.items()},
        samples={t: v[:sample] for t, v in refs.items()},
        pubmed_urls=[
            f"https://pubmed.ncbi.nlm.nih.gov/{p}/" for p in refs["PubMedID"][:sample]
        ],
    )


def pubchem_description(cid: int) -> Description:
    """Short curated descriptions of a compound, with their sources."""
    payload = get_client().compound_json(cid, "description")
    try:
        title, texts = parsers.parse_descriptions(payload)
    except parsers.PayloadError as exc:
        raise PubChemError(f"Unexpected description response: {exc}") from exc
    return Description(cid=cid, title=title, descriptions=texts)


def pubchem_bioactivity(cid: int) -> BioactivitySummary:
    """Summary of PubChem BioAssay results for a CID."""
    try:
        payload = get_client().compound_json(cid, "assaysummary")
    except CompoundNotFound:
        return BioactivitySummary(
            cid=cid,
            tested=False,
            n_assays=0,
            outcomes={},
            top_targets=[],
            most_potent=[],
        )
    try:
        return parsers.summarize_assays(payload, cid)
    except parsers.PayloadError as exc:
        raise PubChemError(f"Unexpected assay summary response: {exc}") from exc


def pubchem_annotations(
    cid: int, heading: AnnotationHeading, max_chars: int = 6000
) -> AnnotationSection:
    """Curated annotations from one section of a compound's PubChem page."""
    if not 500 <= max_chars <= 30000:
        raise BadQuery("max_chars must be between 500 and 30000")
    client = get_client()
    try:
        payload = client.view_data(cid, heading)
    except (CompoundNotFound, BadQuery) as exc:
        try:
            available = parsers.parse_headings(client.view_index(cid))
        except (CompoundNotFound, parsers.PayloadError):
            raise CompoundNotFound(f"No PubChem compound page for CID {cid}") from exc
        shown = ", ".join(available[:40])
        raise AnnotationNotFound(
            f"CID {cid} has no '{heading}' section. Sections on its page include: {shown}",
            available,
        ) from exc
    try:
        title, annotations = parsers.parse_annotations(payload)
    except parsers.PayloadError as exc:
        raise PubChemError(f"Unexpected PUG View response: {exc}") from exc
    kept, truncated = parsers.trim_annotations(annotations, max_chars)
    sources = list(dict.fromkeys(a.source for a in kept if a.source))
    return AnnotationSection(
        cid=cid,
        title=title,
        heading=heading,
        annotations=kept,
        n_total=len(annotations),
        truncated=truncated,
        sources=sources,
    )


def cache_stats(home: Path | str | None = None) -> CacheStats | None:
    """Size of the PubChem response cache, or None if it does not exist."""
    path = cache_path(home)
    return ResponseCache(path).stats() if path.exists() else None


def cache_clear(home: Path | str | None = None) -> int:
    """Empty the PubChem response cache; returns entries removed."""
    path = cache_path(home)
    return ResponseCache(path).clear() if path.exists() else 0


# --------------------------------------------------------------------------- #
# Assistant features
# --------------------------------------------------------------------------- #


def medchem_profile(smiles: str) -> MedchemProfile:
    """Drug-likeness, synthesizability and structural-alert triage."""
    return medchem.profile(parse_smiles(smiles), smiles)


def draw(
    smiles: list[str],
    legends: list[str] | None = None,
    highlight: drawing.Highlight = "none",
    smarts: str | None = None,
    reference: str | None = None,
    size: int = 300,
    per_row: int = 4,
) -> bytes:
    """PNG grid of molecules; see :func:`pubchem_finder.drawing.draw_grid`."""
    return drawing.draw_grid(
        smiles,
        legends=legends,
        highlight=highlight,
        smarts=smarts,
        reference=reference,
        size=size,
        per_row=per_row,
    )
