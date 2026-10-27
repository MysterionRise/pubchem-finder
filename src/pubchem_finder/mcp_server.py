"""MCP server (stdio or streamable HTTP) exposing PubChem lookups and RDKit utilities.

Run with ``pubchem-finder mcp`` (stdio, the default) or
``pubchem-finder mcp --transport http``. Nothing in this module may write to
stdout: over stdio, stdout is the MCP transport.
"""

import functools
import logging
import sys
from collections.abc import Callable
from typing import Literal, TypeVar

import anyio
from mcp.server.mcpserver import Image, MCPServer
from mcp.server.mcpserver.exceptions import ToolError
from starlette.requests import Request
from starlette.responses import JSONResponse, Response

from . import __version__, service
from .core.exceptions import PubChemFinderError
from .data.pug_rest import (
    AnnotationHeading,
    CompoundNotFound,
    IdentityType,
    PubChemTimeout,
    PubChemUnavailable,
)
from .models import (
    AnnotationSection,
    BioactivitySummary,
    Comparison,
    CrossRefs,
    Description,
    IndexStatus,
    LookupResult,
    MedchemProfile,
    Properties,
    RemoteSearchResults,
    SimilarityResults,
    Standardized,
    SubstructureResults,
    Synonyms,
)

T = TypeVar("T")

INSTRUCTIONS = """\
Tools for small-molecule chemistry backed by PubChem and RDKit.
Use pubchem_lookup to resolve names, CIDs, SMILES or InChIKeys to PubChem records
(properties, identifiers, link). Use compound_properties for RDKit descriptors of any
SMILES, standardize to normalize a structure (salts removed, charges neutralized),
and compare for fingerprint similarity between two structures.
similarity_search finds the most similar compounds in the local PubChem index, and
substructure_search the indexed compounds containing a SMILES or SMARTS fragment
(check index_status first; if the index or its substructure search is missing, tell the
user the command it gives).
The pubchem_* tools query PubChem live: structure searches over all of PubChem
(similarity, substructure, exact match), synonyms and CAS numbers, cross-references
(PubMed, patents, registry IDs), descriptions, curated page sections such as GHS
safety and pharmacology, and bioassay summaries. Resolve a name to a CID with
pubchem_lookup first. Cite the sources given in annotation results.
draw_molecules shows structures as an image (highlight a SMARTS pattern, what changed
relative to a reference, or a similarity map); medchem_profile triages drug-likeness,
synthetic accessibility and PAINS/Brenk/NIH structural alerts.
SMILES inputs are case-sensitive; aromatic atoms are lowercase."""

mcp = MCPServer(
    "pubchem-finder",
    version=__version__,
    instructions=INSTRUCTIONS,
    log_level="WARNING",
)

LOOPBACK_HOSTS = frozenset({"127.0.0.1", "localhost", "::1"})


async def _health(_: Request) -> Response:
    """Liveness check for reverse proxies and the HTTP end-to-end test."""
    return JSONResponse(
        {"status": "ok", "name": "pubchem-finder", "version": __version__}
    )


# Registered with the call form: custom_route is untyped as a decorator, which
# mypy --strict rejects. The route only exists over HTTP.
mcp.custom_route("/health", methods=["GET"], include_in_schema=False)(_health)


async def _run(fn: Callable[..., T], *args: object, **kwargs: object) -> T:
    """Run blocking RDKit/HTTP work off the event loop, mapping errors to ToolError."""
    try:
        return await anyio.to_thread.run_sync(functools.partial(fn, *args, **kwargs))
    except CompoundNotFound as exc:
        raise ToolError(
            f"{exc}. Check the spelling or try another identifier."
        ) from exc
    except PubChemUnavailable as exc:
        raise ToolError(f"{exc}. Wait a moment and retry.") from exc
    except PubChemTimeout as exc:
        raise ToolError(
            f"{exc}. Retry in a moment; for structure searches, use a more specific "
            "query, a higher threshold or fewer max_records."
        ) from exc
    except PubChemFinderError as exc:
        raise ToolError(str(exc)) from exc


@mcp.tool()
async def pubchem_lookup(
    query: str,
    by: Literal["name", "cid", "smiles", "inchikey"] = "name",
    limit: int = 5,
) -> LookupResult:
    """Look a compound up in PubChem by name, CID, SMILES or InChIKey.

    Returns PubChem records with CID, title, IUPAC name, SMILES, InChIKey,
    formula, molecular weight, XLogP, H-bond counts, TPSA and the PubChem URL.
    Names resolve to the best matches; use limit to get more than the top 5.
    """
    return await _run(service.lookup, query, by=by, limit=limit)


@mcp.tool()
async def compound_properties(smiles: str) -> Properties:
    """Calculate RDKit descriptors for a SMILES string.

    Returns molecular weight, Crippen LogP, H-bond donors/acceptors, TPSA,
    rotatable bonds, aromatic rings, heavy atoms, canonical SMILES and whether
    Lipinski's Rule of Five passes. Works for any valid SMILES, not only PubChem
    compounds.
    """
    return await _run(service.properties, smiles)


@mcp.tool()
async def standardize(smiles: str) -> Standardized:
    """Standardize a structure: remove salts/counter-ions, neutralize charges,
    normalize functional groups, and return canonical SMILES of the parent."""
    return await _run(service.standardize, smiles)


@mcp.tool()
async def compare(smiles_a: str, smiles_b: str) -> Comparison:
    """Fingerprint similarity between two structures.

    Returns Tanimoto similarity on ECFP4 (Morgan radius 2, 2048 bits) and on
    MACCS keys. 1.0 means identical fingerprints; above ~0.7 on ECFP4 is a close
    analogue.
    """
    return await _run(service.compare, smiles_a, smiles_b)


@mcp.tool()
async def similarity_search(
    smiles: str,
    threshold: float = 0.7,
    top_k: int = 20,
    metric: Literal["tanimoto", "dice", "cosine"] = "tanimoto",
    group_by_scaffold: bool = False,
) -> SimilarityResults:
    """Find the compounds most similar to a structure in the local PubChem index.

    Exact fingerprint similarity (ECFP4/Morgan by default) over every indexed
    compound, best first, with CID, SMILES, InChIKey, IUPAC name, formula, MW
    and PubChem link. threshold is the minimum similarity (0-1, 0.7 finds close
    analogues; lower it to widen the search). Requires a local index built with
    `pubchem-finder index build`; the first search loads it, which can take a
    while for a full PubChem index. group_by_scaffold adds the hits grouped by
    Bemis-Murcko scaffold, to see how many chemical series the hits span.
    """
    return await _run(
        service.similarity_search,
        smiles,
        threshold=threshold,
        top_k=top_k,
        metric=metric,
        group_by_scaffold=group_by_scaffold,
    )


@mcp.tool()
async def substructure_search(
    query: str,
    query_type: Literal["smiles", "smarts"] = "smiles",
    max_hits: int = 50,
) -> SubstructureResults:
    """Find compounds in the local PubChem index that contain a fragment.

    query is a SMILES fragment (e.g. 'c1ccc2[nH]ccc2c1' for indoles) or, with
    query_type="smarts", a SMARTS pattern (e.g. '[NX3][SX4](=O)(=O)c' for aryl
    sulfonamides). Offline, exact (every hit is verified atom by atom; stereo is
    ignored) and not rate-limited, so prefer it over pubchem_substructure_search
    when index_status shows substructure: true. Hits come lowest CID first, which
    puts long-known compounds first. n_candidates bounds how many matches exist;
    exhaustive=false means more matches may exist beyond max_hits (1-1000) or the
    candidate limit. To show where the fragment sits, pass hits to draw_molecules
    with highlight="smarts" and the query as smarts.
    """
    return await _run(
        service.substructure_search, query, query_type=query_type, max_hits=max_hits
    )


@mcp.tool()
async def index_status() -> IndexStatus:
    """Report whether the local PubChem index is built, how many compounds it
    holds, its fingerprint, size, whether substructure search is available, and
    how to build what is missing."""
    return await _run(service.index_status)


@mcp.tool()
async def pubchem_similarity_search(
    smiles: str, threshold: int = 90, max_records: int = 25
) -> RemoteSearchResults:
    """Find similar compounds across all of PubChem (~120M), live.

    Uses PubChem's own 2D fingerprint; threshold is a percentage (0-100, 90 finds
    close analogues). PubChem does not score hits, so each hit gets an ECFP4
    Tanimoto to the query computed here and hits are sorted by it. Prefer the
    local similarity_search when an index is built: it is exact, scored and has
    no rate limits. max_records: 1-200.
    """
    return await _run(
        service.pubchem_similarity_search,
        smiles,
        threshold=threshold,
        max_records=max_records,
    )


@mcp.tool()
async def pubchem_substructure_search(
    smiles: str, max_records: int = 25
) -> RemoteSearchResults:
    """Find compounds across all of PubChem that contain the given structure.

    The query is a SMILES fragment (e.g. 'c1ccc2[nH]ccc2c1' for indoles). Very
    generic fragments can time out; make the query more specific if so. Hits are
    sorted by ECFP4 similarity to the query, closest first. max_records: 1-200.
    """
    return await _run(
        service.pubchem_substructure_search, smiles, max_records=max_records
    )


@mcp.tool()
async def pubchem_exact_match(
    smiles: str, identity_type: IdentityType = "same_stereo_isotope"
) -> RemoteSearchResults:
    """Find the PubChem record(s) for exactly this structure.

    identity_type controls what counts as identical: same_stereo_isotope
    (default, strictest), same_stereo, same_isotope, same_connectivity (ignores
    stereo and isotopes, finds all stereoisomers), same_tautomer,
    nonconflicting_stereo.
    """
    return await _run(service.pubchem_exact_match, smiles, identity_type=identity_type)


@mcp.tool()
async def pubchem_synonyms(cid: int, limit: int = 25) -> Synonyms:
    """Names a compound is known by (common, brand, systematic), most common
    first, plus its CAS Registry Numbers. Get the CID from pubchem_lookup."""
    return await _run(service.pubchem_synonyms, cid, limit=limit)


@mcp.tool()
async def pubchem_xrefs(cid: int) -> CrossRefs:
    """Counts and examples of external identifiers linked to a compound:
    registry IDs (e.g. CAS, ChEMBL, DrugBank), PubMed articles (with links) and
    patents."""
    return await _run(service.pubchem_xrefs, cid)


@mcp.tool()
async def pubchem_description(cid: int) -> Description:
    """Short curated descriptions of a compound (what it is, what it is used
    for), each with its source."""
    return await _run(service.pubchem_description, cid)


@mcp.tool()
async def pubchem_annotations(
    cid: int, heading: AnnotationHeading, max_chars: int = 6000
) -> AnnotationSection:
    """Curated information from one section of a compound's PubChem page.

    Sections: GHS Classification (pictograms, signal word, hazard statements),
    Hazards Identification, Pharmacology and Biochemistry, Mechanism of Action,
    Drug and Medication Information, Experimental Properties (melting point,
    solubility, measured LogP...), Toxicity, Use and Manufacturing, Food Additives
    and Ingredients. Each annotation names its source; cite it. Long sections are
    trimmed to max_chars (500-30000); truncated says when that happened. If the
    page lacks the section, the error lists the sections it has.
    """
    return await _run(service.pubchem_annotations, cid, heading, max_chars=max_chars)


@mcp.tool()
async def pubchem_bioactivity(cid: int) -> BioactivitySummary:
    """Summary of PubChem BioAssay results for a compound: how many assays,
    outcome counts, the targets with the most active results (UniProt
    accession, NCBI gene ID) and the most potent active measurements (µM)."""
    return await _run(service.pubchem_bioactivity, cid)


@mcp.tool()
async def draw_molecules(
    smiles: list[str],
    legends: list[str] | None = None,
    highlight: Literal["none", "smarts", "diff", "similarity_map"] = "none",
    smarts: str | None = None,
    reference_smiles: str | None = None,
) -> Image:
    """Draw 1-24 molecules as a PNG grid image the user can see.

    legends: one caption per molecule (names or CIDs). highlight:
    - "smarts": atoms matching the smarts pattern (e.g. a pharmacophore or an alert),
    - "diff": what each molecule adds relative to reference_smiles (atoms outside
      their maximum common substructure), with drawings aligned to the shared core,
    - "similarity_map": atom contributions to ECFP4 similarity with reference_smiles
      (green raises similarity, pink lowers it).
    Use it to show search hits, analogues or alerts instead of describing structures.
    """
    png = await _run(
        service.draw,
        smiles,
        legends=legends,
        highlight=highlight,
        smarts=smarts,
        reference=reference_smiles,
    )
    return Image(data=png, format="png")


@mcp.tool()
async def medchem_profile(smiles: str) -> MedchemProfile:
    """Medicinal-chemistry triage of a structure: QED drug-likeness (0-1),
    synthetic accessibility (1 easy - 10 hard), Lipinski violations, Veber rules,
    PAINS/Brenk/NIH structural alerts with matched atom indices (pass them to
    draw_molecules with the alert as SMARTS, or just describe them), and the
    Bemis-Murcko scaffold."""
    return await _run(service.medchem_profile, smiles)


def main(transport: str = "stdio", host: str = "127.0.0.1", port: int = 8000) -> None:
    """Serve over stdio (the default) or streamable HTTP at http://HOST:PORT/mcp."""
    logging.basicConfig(
        level=logging.WARNING, format="%(levelname)s %(name)s: %(message)s"
    )
    match transport:
        case "stdio":
            mcp.run(transport="stdio")
        case "http":
            url = f"http://{host}:{port}/mcp"
            # uvicorn announces the address at INFO level, which WARNING hides.
            print(
                f"pubchem-finder {__version__}: MCP over streamable HTTP at {url} "
                "(GET /health to check). Ctrl+C stops it.",
                file=sys.stderr,
                flush=True,
            )
            if host not in LOOPBACK_HOSTS:
                print(
                    "Warning: this server has no authentication or TLS. Keep it on "
                    "127.0.0.1, or put an HTTPS reverse proxy with authentication in "
                    'front of it (docs/clients.md, "Over HTTP").',
                    file=sys.stderr,
                    flush=True,
                )
            mcp.run(transport="streamable-http", host=host, port=port)
        case _:
            raise ValueError(f"Unknown transport {transport!r}; use 'stdio' or 'http'")


if __name__ == "__main__":  # pragma: no cover
    main()
