"""MCP server (stdio) exposing PubChem lookups and RDKit utilities.

Run with ``pubchem-finder mcp``. Nothing in this module may write to stdout:
stdout is the MCP transport.
"""

import functools
import logging
from collections.abc import Callable
from typing import Literal, TypeVar

import anyio
from mcp.server.mcpserver import MCPServer
from mcp.server.mcpserver.exceptions import ToolError

from . import __version__, service
from .core.exceptions import PubChemFinderError
from .data.pug_rest import CompoundNotFound, PubChemUnavailable
from .models import (
    Comparison,
    IndexStatus,
    LookupResult,
    Properties,
    SimilarityResults,
    Standardized,
)

T = TypeVar("T")

INSTRUCTIONS = """\
Tools for small-molecule chemistry backed by PubChem and RDKit.
Use pubchem_lookup to resolve names, CIDs, SMILES or InChIKeys to PubChem records
(properties, identifiers, link). Use compound_properties for RDKit descriptors of any
SMILES, standardize to normalize a structure (salts removed, charges neutralized),
and compare for fingerprint similarity between two structures.
similarity_search finds the most similar compounds in the local PubChem index
(check index_status first; if the index is not built, tell the user how to build it).
SMILES inputs are case-sensitive; aromatic atoms are lowercase."""

mcp = MCPServer(
    "pubchem-finder",
    version=__version__,
    instructions=INSTRUCTIONS,
    log_level="WARNING",
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
) -> SimilarityResults:
    """Find the compounds most similar to a structure in the local PubChem index.

    Exact fingerprint similarity (ECFP4/Morgan by default) over every indexed
    compound, best first, with CID, SMILES, InChIKey, IUPAC name, formula, MW
    and PubChem link. threshold is the minimum similarity (0-1, 0.7 finds close
    analogues; lower it to widen the search). Requires a local index built with
    `pubchem-finder index build`; the first search loads it, which can take a
    while for a full PubChem index.
    """
    return await _run(
        service.similarity_search,
        smiles,
        threshold=threshold,
        top_k=top_k,
        metric=metric,
    )


@mcp.tool()
async def index_status() -> IndexStatus:
    """Report whether the local PubChem similarity index is built, how many
    compounds it holds, its fingerprint, size, and how to build it if missing."""
    return await _run(service.index_status)


def main() -> None:
    """Serve over stdio."""
    logging.basicConfig(
        level=logging.WARNING, format="%(levelname)s %(name)s: %(message)s"
    )
    mcp.run(transport="stdio")


if __name__ == "__main__":  # pragma: no cover
    main()
