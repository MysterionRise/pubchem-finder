"""MCP server (stdio or streamable HTTP) exposing PubChem lookups and RDKit utilities.

Run with ``pubchem-finder mcp`` (stdio, the default) or
``pubchem-finder mcp --transport http``. Nothing in this module may write to
stdout: over stdio, stdout is the MCP transport.
"""

import functools
import logging
import sys
from collections.abc import Callable
from typing import Any, Literal, TypeVar

import anyio
from mcp.server.mcpserver import Image, MCPServer
from mcp.server.mcpserver.exceptions import (
    ResourceError,
    ResourceNotFoundError,
    ToolError,
)
from mcp.shared.exceptions import MCPError
from mcp.types import INVALID_PARAMS, Icon, ToolAnnotations
from pydantic import BaseModel
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
from .data.sdf import parse_id
from .models import (
    AnnotationSection,
    BioactivitySummary,
    Comparison,
    CompoundRecord,
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

logger = logging.getLogger(__name__)

T = TypeVar("T")
F = TypeVar("F", bound=Callable[..., Any])

INSTRUCTIONS = """\
Tools for small-molecule chemistry backed by PubChem and RDKit.
Use pubchem_lookup to resolve names, CIDs, SMILES or InChIKeys to PubChem records
(properties, identifiers, link). Use compound_properties for RDKit descriptors of any
SMILES, standardize to normalize a structure (salts removed, charges neutralized),
and compare for fingerprint similarity between two structures.
similarity_search finds the most similar compounds in the local index, and
substructure_search the indexed compounds containing a SMILES or SMARTS fragment.
Until the user builds a local PubChem index, both search a built-in demo set of about
1,900 FDA-approved drugs and say so in note: tell the user, and offer the pubchem_*
searches (all of PubChem, live) or the build command index_status gives.
The pubchem_* tools query PubChem live: structure searches over all of PubChem
(similarity, substructure, exact match), synonyms and CAS numbers, cross-references
(PubMed, patents, registry IDs), descriptions, curated page sections such as GHS
safety and pharmacology, and bioassay summaries. Resolve a name to a CID with
pubchem_lookup first. Cite the sources given in annotation results.
draw_molecules shows structures as an image (highlight a SMARTS pattern, what changed
relative to a reference, or a similarity map); medchem_profile triages drug-likeness,
synthetic accessibility and PAINS/Brenk/NIH structural alerts.
SMILES inputs are case-sensitive; aromatic atoms are lowercase."""

REPO_URL = "https://github.com/MysterionRise/pubchem-finder"
# The plugin's icon, served from the repository rather than inlined: a data: URI
# would add ~35 KB to every initialize response.
ICON_URL = (
    "https://raw.githubusercontent.com/MysterionRise/pubchem-finder/master/"
    "plugins/pubchem-finder/.claude-plugin/icon.png"
)

mcp = MCPServer(
    "pubchem-finder",
    title="PubChem Finder",
    version=__version__,
    instructions=INSTRUCTIONS,
    website_url=REPO_URL,
    icons=[Icon(src=ICON_URL, mime_type="image/png", sizes=["1024x1024"])],
    log_level="WARNING",
)


# Claude Code saves a text result longer than 50,000 characters to a file instead of
# showing it, unless the tool raises its own limit (up to 500,000). Searches with many
# hits and long annotation sections can pass 50,000; this lets them through whole.
LARGE_RESULTS = {"anthropic/maxResultSizeChars": 100_000}


def _tool(
    title: str, *, calls_pubchem: bool, large_results: bool = False
) -> Callable[[F], F]:
    """Register a tool with a title and annotations, as directory policies require.

    Every tool only reads: its arguments, the local index or PubChem. The response
    cache and the demo index (built on the first local search when no index
    exists) are internal state, not the tool's environment. ``calls_pubchem`` marks
    the tools that reach an external service (the open world); ``large_results``
    the ones whose results can pass Claude Code's default size limit.
    """
    return mcp.tool(
        title=title,
        annotations=ToolAnnotations(
            title=title,
            read_only_hint=True,
            destructive_hint=False,
            idempotent_hint=True,
            open_world_hint=calls_pubchem,
        ),
        meta=dict(LARGE_RESULTS) if large_results else None,
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


# PubChem structure searches that time out because the query matches too much; for
# these, narrowing the query helps, while for lookups only a retry does.
_NARROWABLE_SEARCHES = frozenset(
    {"pubchem_similarity_search", "pubchem_substructure_search"}
)


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
        hint = "Retry in a moment"
        if getattr(fn, "__name__", "") in _NARROWABLE_SEARCHES:
            hint += (
                "; use a more specific query, a higher threshold or fewer max_records"
            )
        raise ToolError(f"{exc}. {hint}.") from exc
    except PubChemFinderError as exc:
        raise ToolError(str(exc)) from exc
    except Exception as exc:  # a damaged cache or index file, a disk error, or a bug
        logger.exception("%s failed", getattr(fn, "__name__", "tool"))
        detail = f"{type(exc).__name__}: {exc}"[:300]
        raise ToolError(
            f"Unexpected error ({detail}). If it persists, report it at "
            f"{REPO_URL}/issues with the server's log."
        ) from exc


@_tool("Look up a PubChem compound", calls_pubchem=True)
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


@_tool("Compute molecular properties", calls_pubchem=False)
async def compound_properties(smiles: str) -> Properties:
    """Calculate RDKit descriptors for a SMILES string.

    Returns molecular weight, Crippen LogP, H-bond donors/acceptors, TPSA,
    rotatable bonds, aromatic rings, heavy atoms, canonical SMILES and whether
    Lipinski's Rule of Five passes. Works for any valid SMILES, not only PubChem
    compounds.
    """
    return await _run(service.properties, smiles)


@_tool("Standardize a structure", calls_pubchem=False)
async def standardize(smiles: str) -> Standardized:
    """Standardize a structure: remove salts/counter-ions, neutralize charges,
    normalize functional groups, and return canonical SMILES of the parent.

    The parent is the largest organic component; `removed` lists everything left
    out, so for a mixture of two drugs it names the one that was dropped. Metal
    compounds (cisplatin) and inputs with nothing organic (sodium chloride) are
    returned as given, with a `note`. Sulfoxides stay S=O, the form PubChem and
    the local index use."""
    return await _run(service.standardize, smiles)


@_tool("Compare two structures", calls_pubchem=False)
async def compare(smiles_a: str, smiles_b: str) -> Comparison:
    """Fingerprint similarity between two structures.

    Returns Tanimoto similarity on ECFP4 (Morgan radius 2, 2048 bits) and on
    MACCS keys. 1.0 means identical fingerprints; above ~0.7 on ECFP4 is a close
    analogue.
    """
    return await _run(service.compare, smiles_a, smiles_b)


@_tool("Similarity search (local index)", calls_pubchem=False, large_results=True)
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
    analogues; lower it to widen the search). Searches the index built with
    `pubchem-finder index build` (the first search loads it, which takes a while
    for a full PubChem index); until one is built, a demo set of about 1,900
    FDA-approved drugs, with search_scope="demo" and a note to pass on to the
    user. group_by_scaffold adds the hits grouped by
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


@_tool("Substructure search (local index)", calls_pubchem=False, large_results=True)
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
    when index_status shows substructure: true. Without a local index it
    searches the demo set of FDA-approved drugs (search_scope="demo"). Hits come
    lowest CID first, which puts long-known compounds first. n_candidates bounds
    how many matches exist; exhaustive=false means more matches may exist beyond
    max_hits (1-1000) or the candidate limit. To show where the fragment sits,
    pass hits to draw_molecules with highlight="smarts" and the query as smarts.
    """
    return await _run(
        service.substructure_search, query, query_type=query_type, max_hits=max_hits
    )


@_tool("Local index status", calls_pubchem=False)
async def index_status() -> IndexStatus:
    """Report whether the local PubChem index is built, how many compounds it
    holds, its fingerprint, size, whether substructure search is available, what
    local searches cover now (search_scope: the index, or the built-in demo set of
    FDA-approved drugs) and how to build what is missing."""
    return await _run(service.index_status)


@_tool("PubChem similarity search", calls_pubchem=True)
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


@_tool("PubChem substructure search", calls_pubchem=True)
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


@_tool("PubChem exact structure match", calls_pubchem=True)
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


@_tool("PubChem synonyms and CAS numbers", calls_pubchem=True)
async def pubchem_synonyms(cid: int, limit: int = 25) -> Synonyms:
    """Names a compound is known by (common, brand, systematic), most common
    first, plus its CAS Registry Numbers. Get the CID from pubchem_lookup."""
    return await _run(service.pubchem_synonyms, cid, limit=limit)


@_tool("PubChem cross-references", calls_pubchem=True)
async def pubchem_xrefs(cid: int) -> CrossRefs:
    """Counts and examples of external identifiers linked to a compound:
    registry IDs (e.g. CAS, ChEMBL, DrugBank), PubMed articles (with links) and
    patents."""
    return await _run(service.pubchem_xrefs, cid)


@_tool("PubChem compound description", calls_pubchem=True)
async def pubchem_description(cid: int) -> Description:
    """Short curated descriptions of a compound (what it is, what it is used
    for), each with its source."""
    return await _run(service.pubchem_description, cid)


@_tool("PubChem page section", calls_pubchem=True, large_results=True)
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


@_tool("PubChem bioactivity summary", calls_pubchem=True)
async def pubchem_bioactivity(cid: int) -> BioactivitySummary:
    """Summary of PubChem BioAssay results for a compound: how many assays,
    outcome counts, the targets with the most active assays (one per gene, with
    its NCBI gene ID and a protein accession) and the lowest potency-type values
    among active results (IC50, EC50, AC50, Ki, Kd and similar, in µM). The lowest
    value is not always the most meaningful one: name the assay when you quote it."""
    return await _run(service.pubchem_bioactivity, cid)


@_tool("Draw molecules", calls_pubchem=False)
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


@_tool("Medicinal chemistry profile", calls_pubchem=False)
async def medchem_profile(smiles: str) -> MedchemProfile:
    """Medicinal-chemistry triage of a structure: QED drug-likeness (0-1),
    synthetic accessibility (1 easy - 10 hard), Lipinski violations, Veber rules,
    PAINS/Brenk/NIH structural alerts with matched atom indices (pass them to
    draw_molecules with the alert as SMARTS, or just describe them), and the
    Bemis-Murcko scaffold."""
    return await _run(service.medchem_profile, smiles)


# --------------------------------------------------------------------------- #
# Resources: data a client can attach to a conversation (in Claude Code,
# @pubchem-finder:pubchem://compound/2244).
# --------------------------------------------------------------------------- #


async def _read(fn: Callable[..., BaseModel], *args: object) -> str:
    """Run a service call off the event loop and return its result as JSON; every
    failure becomes a ResourceError with the reason (the SDK's default hides it)."""
    try:
        result = await anyio.to_thread.run_sync(functools.partial(fn, *args))
    except CompoundNotFound as exc:
        raise ResourceNotFoundError(str(exc)) from exc
    except PubChemFinderError as exc:
        raise ResourceError(str(exc)) from exc
    except Exception as exc:  # see _run
        logger.exception("reading a resource failed")
        raise ResourceError(
            f"Unexpected error ({type(exc).__name__}: {exc})"[:300]
        ) from exc
    return result.model_dump_json(indent=2)


@mcp.resource(
    "pubchem://compound/{cid}",
    title="PubChem compound",
    description="A PubChem compound record by CID: names, SMILES, InChIKey, formula, "
    "weight, XLogP, H-bond counts, TPSA and link (fetched live, then cached).",
    mime_type="application/json",
)
async def compound_resource(cid: str) -> str:
    number = parse_id(cid)  # ASCII digits, leading zeros allowed, at most 32 bits
    if number is None:
        raise ResourceNotFoundError(f"Not a PubChem CID: {cid[:20]!r}")
    return await _read(_compound, number)


def _compound(cid: int) -> CompoundRecord:
    return service.lookup(str(cid), by="cid", limit=1).compounds[0]


@mcp.resource(
    "pubchem://index/status",
    title="Local index status",
    description="Whether the local index is built, what it holds and what local "
    "searches cover now (the index, or the built-in demo set).",
    mime_type="application/json",
)
async def index_status_resource() -> str:
    return await _read(service.index_status)


# --------------------------------------------------------------------------- #
# Prompts: workflows people start by name (in Claude Code,
# /pubchem-finder:compound_dossier aspirin). Claude Code splits prompt arguments
# on whitespace, so each one is a single token.
# --------------------------------------------------------------------------- #


@mcp.prompt(title="Compound dossier")
def compound_dossier(compound: str) -> str:
    """A cited profile of one compound: identity, structure, drug-likeness, safety,
    pharmacology and bioassay targets. compound: a name, CID or SMILES."""
    return (
        f"Write a concise dossier on {compound} from PubChem, with a source for every "
        "fact.\n"
        '1. Resolve it to one PubChem compound with pubchem_lookup (by="cid" for a '
        'number, by="smiles" for a structure); if several match, say which you chose.\n'
        "2. Identity: pubchem_synonyms (common names, CAS numbers) and "
        "pubchem_description.\n"
        "3. Structure: draw_molecules, compound_properties and medchem_profile on its "
        "SMILES.\n"
        '4. Safety: pubchem_annotations with heading "GHS Classification".\n'
        '5. Pharmacology: pubchem_annotations with "Pharmacology and Biochemistry" and '
        '"Mechanism of Action", and pubchem_bioactivity for its assay targets.\n'
        "When its PubChem page lacks a section, say so instead of filling it in from "
        "memory. End with the PubChem link."
    )


@mcp.prompt(title="Find analogues")
def find_analogues(compound: str, threshold: str = "0.6") -> str:
    """Close structural analogues of a compound, grouped into series and drawn with
    what changed. compound: a name, CID or SMILES; threshold: the minimum ECFP4
    Tanimoto similarity, 0-1 (default 0.6)."""
    try:
        value = float(threshold)
    except ValueError:
        value = -1.0
    if not 0.0 < value <= 1.0:
        raise MCPError(
            INVALID_PARAMS, f"threshold must be a number in (0, 1], not {threshold!r}"
        )
    return (
        f"Find close analogues of {compound}.\n"
        "1. Get its SMILES with pubchem_lookup unless it is one already.\n"
        f"2. Run similarity_search with threshold={value:g} and group_by_scaffold=true. "
        'If the result\'s search_scope is "demo", tell me that only the built-in demo '
        "set of FDA-approved drugs was searched, and also run pubchem_similarity_search "
        "to cover all of PubChem.\n"
        "3. Draw the closest hits (up to 8, the query first) with draw_molecules, "
        'highlight="diff" and the query as reference_smiles, with CIDs and similarities '
        "as legends.\n"
        "4. Summarize how many series (scaffolds) the hits span, what the closest "
        "analogues change, and which are known drugs."
    )


@mcp.prompt(title="Triage structures")
def triage_structures(structures: str = "") -> str:
    """Medicinal-chemistry triage of a set of structures, ranked: drug-likeness,
    synthetic accessibility and structural alerts. structures: SMILES separated by
    commas without spaces, or empty to give them afterwards."""
    items = [s.strip() for s in structures.split(",") if s.strip()]
    target = (
        "these structures:\n" + "\n".join(f"- {s}" for s in items)
        if items
        else "the structures I give you next (ask me for them as SMILES)."
    )
    return (
        f"Triage {target}\n"
        "For each one, standardize it (salts, charges; say what it removed when that "
        "is more than a counter-ion, and pass its note on), then run medchem_profile: QED, "
        "synthetic accessibility, Lipinski and Veber rules, PAINS/Brenk/NIH alerts. Show "
        "one table ranked from most to least promising, with the reason for each "
        "rank. Then draw the structures that have alerts with draw_molecules, their "
        "alert names as legends, and say what each alert is about."
    )


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
