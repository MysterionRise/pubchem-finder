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

from . import drawing, hints, medchem
from .chemistry import MoleculeStandardizer, calculate_properties, is_druglike
from .config import cache_enabled, cache_path, demo_enabled, resolve_home
from .core.exceptions import MoleculeValidationError, PubChemFinderError
from .core.limits import check_atoms, check_length
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
    PubChemTimeout,
    PubChemUnavailable,
)
from .fingerprints import MACCSFingerprint, MorganFingerprint
from .index import demo
from .index.builder import DEFAULT_FP_PARAMS
from .index.engine import LocalIndex, Metric
from .index.layout import Merged
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
# PubChem's hits are re-scored with the local index's default fingerprint, so
# their scores are comparable with similarity_search; compare() keeps 2048 bits.
_ECFP4_INDEX = MorganFingerprint(radius=2, n_bits=int(DEFAULT_FP_PARAMS["fpSize"]))
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
    """Parse SMILES or raise ``MoleculeValidationError``, also for oversized input
    (see :mod:`pubchem_finder.core.limits`)."""
    check_length(smiles.strip())
    mol: Chem.Mol | None = (
        Chem.MolFromSmiles(smiles.strip()) if smiles.strip() else None
    )
    if mol is None:
        raise MoleculeValidationError(f"Invalid SMILES: {smiles!r}")
    return check_atoms(mol)


def parse_query(query: str, query_type: QueryType = "smiles") -> Chem.Mol:
    """Parse a substructure query given as SMILES or SMARTS."""
    if query_type == "smiles":
        return parse_smiles(query)
    if query_type != "smarts":
        raise PubChemFinderError(f"Unknown query_type {query_type!r}")
    check_length(query.strip(), "SMARTS")
    mol: Chem.Mol | None = Chem.MolFromSmarts(query.strip()) if query.strip() else None
    if mol is None:
        raise MoleculeValidationError(f"Invalid SMARTS: {query!r}")
    check_atoms(mol, "SMARTS")
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
    std = _STANDARDIZER.standardize_detailed(mol)
    if std is None:
        raise MoleculeValidationError(f"Could not standardize {smiles!r}")
    result = str(Chem.MolToSmiles(std.mol))
    return Standardized(
        input_smiles=smiles,
        smiles=result,
        changed=result != Chem.MolToSmiles(mol),
        removed=std.removed,
        note=std.note,
    )


def compare(smiles_a: str, smiles_b: str) -> Comparison:
    """Tanimoto similarity between two structures on ECFP4 (Morgan radius 2,
    2048 bits) and on MACCS keys."""
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


DEMO_NOTE = (
    "Searched the built-in demo set of {n:,} FDA-approved drugs, not PubChem. "
    "pubchem_similarity_search and pubchem_substructure_search (`{live}` in a "
    "terminal) search all of PubChem live. No local index is built yet. {build}"
)


def _search_index(
    home: Path | str | None, live: str = "pubchem similarity SMILES"
) -> tuple[LocalIndex, str, str | None]:
    """The index a local search uses, its scope and the note for the results.
    ``live`` is the command that runs the same search over all of PubChem.

    The local index once one is built. Until then the built-in demo index (built
    on first use, in about a second), unless PUBCHEM_FINDER_DEMO=0; then the
    local index anyway, whose search says how to build it.
    """
    index = get_index(home)
    if not demo_enabled() or index.paths.load_manifest() is not None:
        return index, "index", None
    demo_index = get_index(demo.ensure_built(index.paths.home))
    note = DEMO_NOTE.format(
        n=demo_index.manifest.n_compounds,
        live=hints.command(live),
        build=hints.build_index(index.paths.home),
    )
    return demo_index, "demo", note


_UNMERGED_EFFECT = {
    "similarity": "hits from there come without a structure",
    "substructure": "its compounds cannot be checked and are left out",
    "status": "similarity hits from there come without a structure, substructure "
    "search leaves its compounds out",
}


def _unmerged_warning(
    index: LocalIndex,
    merged: Merged | None,
    kind: Literal["similarity", "substructure", "status"],
) -> str | None:
    """What to tell the user when chunks changed after the index was merged.

    ``merged`` is the merge the answer came from (None: the index on disk now).
    """
    stale = index.unmerged(merged)
    if not stale:
        return None

    def some(names: list[str]) -> str:
        more = f" and {len(names) - 3} more" if len(names) > 3 else ""
        return ", ".join(names[:3]) + more

    rebuilt = [name for name, stats in stale.items() if stats is not None]
    missing = [name for name, stats in stale.items() if stats is None]
    what = []
    if rebuilt:
        what.append(
            f"{some(rebuilt)} was indexed again after the last merge (a build "
            "failed, was stopped, or is still running)"
        )
    if missing:
        what.append(
            f"{some(missing)} is no longer complete in the chunks directory (a "
            "build of it is running or was stopped, or its files were removed)"
        )
    files = [
        stats.source.path if stats is not None and stats.source is not None else None
        for stats in stale.values()
    ]
    return (
        f"The index is behind its files: {'; '.join(what)}. Until a build "
        f"completes, {_UNMERGED_EFFECT[kind]}, and compounds added since are not "
        f"found. Finish with {hints.merge_rebuilt(files, index.paths.home)}."
    )


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
    index, scope, note = _search_index(home)
    # The hits are looked up, and judged, against the merge that gave them: a
    # build can finish between the two steps.
    pairs, merged = index.search_merged(
        mol, threshold=threshold, top_k=top_k, metric=metric
    )
    manifest = merged.manifest
    meta = index.metadata([cid for cid, _ in pairs], merged)
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
        search_scope=scope,
        note=note,
        hits=hits,
        scaffolds=medchem.group_by_scaffold(hits) if group_by_scaffold else None,
        warning=_unmerged_warning(index, merged, "similarity"),
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
    index, scope, note = _search_index(home, live="pubchem substructure SMILES")
    found = index.substructure(mol, max_hits=max_hits, max_candidates=max_candidates)
    sub = found.merged.manifest.substructure
    warnings = [_unmerged_warning(index, found.merged, "substructure")]
    if sub is not None and sub.rdkit_version != rdkit.__version__:
        warnings.append(
            f"The substructure index was built with RDKit {sub.rdkit_version}, but "
            f"RDKit {rdkit.__version__} is running. Pattern fingerprints can differ "
            "between releases, so matches may be missing. Rebuild it with "
            f"{hints.add_substructure(index.paths.home)}."
        )
    warning = " ".join(w for w in warnings if w) or None
    return SubstructureResults(
        query=query,
        query_type=query_type,
        n_indexed=sub.n_compounds if sub is not None else 0,
        search_scope=scope,
        note=note,
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
    merged = paths.read_manifest()
    manifest = merged.manifest if merged is not None else None
    demo_compounds = int(demo.provenance()["n_compounds"]) if demo_enabled() else 0
    if manifest is None:
        until_then = (
            f" Until then, similarity_search and substructure_search search the "
            f"built-in demo set of {demo_compounds:,} FDA-approved drugs."
            if demo_compounds
            else ""
        )
        return IndexStatus(
            state="not_built",
            home=str(paths.home),
            chunks_built=len(built),
            search_scope="demo" if demo_compounds else "none",
            demo_compounds=demo_compounds,
            hint=hints.build_index(paths.home) + until_then,
        )
    sub = manifest.substructure if paths.pattern_fingerprints.exists() else None
    hint_parts = [
        _unmerged_warning(index, merged, "status"),
        None
        if sub is not None
        else (
            f"Substructure search is off; add it with {hints.add_substructure(paths.home)}."
        ),
    ]
    if manifest.chunk_builds is None and any(c.source is not None for c in built):
        hint_parts.insert(
            1,
            "This index was merged by pubchem-finder 2.2.1 or earlier, which could "
            "miss a file that had been indexed again; the next `index build` merges "
            "it once more to be sure.",
        )
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
        unmerged_chunks=list(index.unmerged(merged)),
        loaded=index.loaded,
        mode=None
        if index.in_memory is None
        else ("memory" if index.in_memory else "disk"),
        substructure=sub is not None,
        substructure_size_bytes=paths.pattern_fingerprints.stat().st_size
        if sub is not None
        else 0,
        substructure_memory_bytes=sub.memory_bytes() if sub is not None else 0,
        search_scope="index",
        demo_compounds=demo_compounds,
        hint=" ".join(part for part in hint_parts if part) or None,
    )


# --------------------------------------------------------------------------- #
# PubChem knowledge tools
# --------------------------------------------------------------------------- #

SIMILARITY_NOTE = (
    "PubChem selected these compounds with its own 881-bit substructure fingerprint "
    "(Tanimoto >= {threshold}%), which does not report a score per hit. ecfp4_similarity "
    "is computed here on ECFP4 (Morgan radius 2, 1024 bits, the local index default), "
    "so it is comparable with similarity_search; hits are sorted by it."
)
SUBSTRUCTURE_NOTE = (
    "Every hit contains the query as a substructure. PubChem returns at most max_records "
    "hits; they are sorted by ecfp4_similarity to the query, so close analogues come first."
)
IDENTITY_NOTE = "Compounds PubChem considers identical to the query under identity_type={identity_type}."


def _rank(query: Chem.Mol, records: list[CompoundRecord]) -> list[RemoteHit]:
    query_fp = _ECFP4_INDEX.generate_bitvect(query)
    hits: list[RemoteHit] = []
    for record in records:
        similarity = None
        smiles = record.smiles or record.connectivity_smiles
        mol = Chem.MolFromSmiles(smiles) if smiles else None
        if mol is not None:
            similarity = round(
                float(
                    DataStructs.TanimotoSimilarity(
                        query_fp, _ECFP4_INDEX.generate_bitvect(mol)
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


# Patents get a request of their own: for well-studied compounds they are most of the
# response (3.1 of 3.6 MB for aspirin) and the part PubChem runs out of time on.
_XREF_PATENTS = ("PatentID",)
_XREF_OTHERS = tuple(t for t in XREF_TYPES if t not in _XREF_PATENTS)


def _require_compound(cid: int) -> None:
    """Raise ``CompoundNotFound`` unless PubChem has ``cid``.

    PubChem answers "no such compound" and "a compound without this data" with the
    same 404, so an empty result is checked once; the answer is cached.
    """
    get_client().lookup(str(cid), by="cid", limit=1)


def _xrefs(cid: int, types: tuple[str, ...]) -> dict[str, list[str]]:
    try:
        payload = get_client().compound_json(cid, "xrefs/" + ",".join(types))
        return parsers.parse_xrefs(payload, types)
    except CompoundNotFound:  # PubChem's answer when a compound has none of these
        return {t: [] for t in types}
    except parsers.PayloadError as exc:
        raise PubChemError(f"Unexpected xrefs response: {exc}") from exc


def pubchem_xrefs(cid: int, sample: int = 10) -> CrossRefs:
    """Registry IDs, PubMed articles and patents linked to a CID.

    When only the patent request fails (PubChem runs out of time, is busy, or its
    answer cannot be read), the other types are still returned, patents are listed
    in ``unavailable`` and ``note`` says why.
    """
    if not 1 <= sample <= 100:
        raise BadQuery("sample must be between 1 and 100")
    refs = _xrefs(cid, _XREF_OTHERS)
    unavailable: list[str] = []
    note: str | None = None
    try:
        refs |= _xrefs(cid, _XREF_PATENTS)
    except BadQuery:
        raise
    except PubChemError as exc:
        logger.warning("patents for CID %d left out: %s", cid, exc)
        unavailable = list(_XREF_PATENTS)
        note = _patents_note(exc)
    if not any(refs.values()):
        _require_compound(cid)
    return CrossRefs(
        cid=cid,
        counts={t: len(v) for t, v in refs.items()},
        samples={t: v[:sample] for t, v in refs.items()},
        pubmed_urls=[
            f"https://pubmed.ncbi.nlm.nih.gov/{p}/" for p in refs["PubMedID"][:sample]
        ],
        unavailable=unavailable,
        note=note,
    )


def _patents_note(exc: PubChemError) -> str:
    rest = "The other identifiers are complete."
    if isinstance(exc, PubChemTimeout):
        return (
            "PubChem ran out of time listing the patents. It stops after 30 seconds, "
            "and for the most-patented compounds (water, ethanol) the list takes "
            f"longer every time, so asking again may not help. {rest}"
        )
    if isinstance(exc, PubChemUnavailable):
        return f"PubChem was too busy to list the patents; ask again later. {rest}"
    return f"The patent list could not be read ({exc}); ask again later. {rest}"


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
        _require_compound(cid)  # untested, or not a compound at all?
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
