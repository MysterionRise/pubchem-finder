"""Command-line interface: ``pubchem-finder``."""

import sys
from enum import StrEnum
from typing import Annotated

import typer
from pydantic import BaseModel

from . import __version__, service
from .core.exceptions import PubChemFinderError

app = typer.Typer(
    help="Molecular similarity search over PubChem: lookups, descriptors, MCP server.",
    no_args_is_help=True,
    rich_markup_mode=None,
)


class By(StrEnum):
    name = "name"
    cid = "cid"
    smiles = "smiles"
    inchikey = "inchikey"


JsonFlag = Annotated[bool, typer.Option("--json", help="Emit JSON instead of text.")]


def _version(value: bool) -> None:
    if value:
        typer.echo(f"pubchem-finder {__version__}")
        raise typer.Exit()


@app.callback()
def _main(
    version: Annotated[
        bool,
        typer.Option(
            "--version", callback=_version, is_eager=True, help="Show version."
        ),
    ] = False,
) -> None:
    pass


def _emit(model: BaseModel, as_json: bool) -> None:
    if as_json:
        typer.echo(model.model_dump_json(indent=2))
        return
    for key, value in model.model_dump().items():
        if isinstance(value, list):
            continue
        if isinstance(value, float):
            value = f"{value:.3f}"
        typer.echo(f"{key:<20} {value}")


def _fail(exc: Exception) -> None:
    typer.echo(f"error: {exc}", err=True)
    raise typer.Exit(code=1)


@app.command()
def lookup(
    query: Annotated[str, typer.Argument(help="Name, CID, SMILES or InChIKey.")],
    by: Annotated[By, typer.Option(help="Identifier type of QUERY.")] = By.name,
    limit: Annotated[int, typer.Option(min=1, help="Max records.")] = 5,
    as_json: JsonFlag = False,
) -> None:
    """Look a compound up in PubChem."""
    try:
        result = service.lookup(query, by=by.value, limit=limit)
    except PubChemFinderError as exc:
        _fail(exc)
        return
    if as_json:
        typer.echo(result.model_dump_json(indent=2))
        return
    for record in result.compounds:
        typer.echo(f"CID {record.cid}  {record.title or ''}")
        typer.echo(f"  SMILES   {record.smiles or ''}")
        typer.echo(
            f"  Formula  {record.molecular_formula or ''}  MW {record.molecular_weight}"
        )
        typer.echo(f"  InChIKey {record.inchikey or ''}")
        typer.echo(f"  {record.pubchem_url}")


@app.command()
def props(
    smiles: Annotated[str, typer.Argument(help="SMILES string.")],
    as_json: JsonFlag = False,
) -> None:
    """RDKit descriptors and Lipinski check for a structure."""
    try:
        _emit(service.properties(smiles), as_json)
    except PubChemFinderError as exc:
        _fail(exc)


@app.command()
def standardize(
    smiles: Annotated[str, typer.Argument(help="SMILES string.")],
    as_json: JsonFlag = False,
) -> None:
    """Standardize a structure (remove salts, neutralize, canonicalize)."""
    try:
        result = service.standardize(smiles)
    except PubChemFinderError as exc:
        _fail(exc)
        return
    if as_json:
        typer.echo(result.model_dump_json(indent=2))
    else:
        typer.echo(result.smiles)


@app.command()
def compare(
    smiles_a: Annotated[str, typer.Argument(help="First SMILES.")],
    smiles_b: Annotated[str, typer.Argument(help="Second SMILES.")],
    as_json: JsonFlag = False,
) -> None:
    """Fingerprint similarity (ECFP4 and MACCS Tanimoto) between two structures."""
    try:
        _emit(service.compare(smiles_a, smiles_b), as_json)
    except PubChemFinderError as exc:
        _fail(exc)


@app.command()
def mcp() -> None:
    """Run the MCP server over stdio (for Claude Desktop, Claude Code, etc.)."""
    from .mcp_server import main

    main()


def run() -> None:  # pragma: no cover - console script entry
    sys.exit(app())


if __name__ == "__main__":  # pragma: no cover
    run()
