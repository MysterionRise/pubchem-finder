"""Command-line interface: ``pubchem-finder``."""

import os
import sys
from enum import StrEnum
from pathlib import Path
from typing import Annotated

import typer
from pydantic import BaseModel

from . import __version__, service
from .config import resolve_home
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


class Metric(StrEnum):
    tanimoto = "tanimoto"
    dice = "dice"
    cosine = "cosine"


JsonFlag = Annotated[bool, typer.Option("--json", help="Emit JSON instead of text.")]
HomeOpt = Annotated[
    Path | None,
    typer.Option(
        "--home",
        envvar="PUBCHEM_FINDER_HOME",
        help="Data directory (default: platform user data dir).",
        show_default=False,
    ),
]

index_app = typer.Typer(
    help="Build and inspect the local similarity index.",
    no_args_is_help=True,
    rich_markup_mode=None,
)
app.add_typer(index_app, name="index")


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
        if value is None or isinstance(value, list):
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
def search(
    smiles: Annotated[str, typer.Argument(help="Query SMILES.")],
    threshold: Annotated[
        float, typer.Option(min=0.01, max=1.0, help="Minimum similarity (0-1).")
    ] = 0.7,
    top: Annotated[int, typer.Option(min=1, max=1000, help="Max hits.")] = 20,
    metric: Annotated[
        Metric, typer.Option(help="Similarity metric.")
    ] = Metric.tanimoto,
    home: HomeOpt = None,
    as_json: JsonFlag = False,
) -> None:
    """Find the most similar compounds in the local index."""
    try:
        result = service.similarity_search(
            smiles, threshold=threshold, top_k=top, metric=metric.value, home=home
        )
    except PubChemFinderError as exc:
        _fail(exc)
        return
    if as_json:
        typer.echo(result.model_dump_json(indent=2))
        return
    typer.echo(
        f"{len(result.hits)} hit(s) >= {result.threshold} {result.metric} "
        f"among {result.n_indexed:,} compounds ({result.fingerprint})"
    )
    for hit in result.hits:
        name = hit.iupac_name or ""
        typer.echo(
            f"{hit.similarity:.3f}  CID {hit.cid:<10} {hit.smiles or ''}  {name}"
        )


def _default_workers() -> int:
    return max(1, (os.cpu_count() or 2) - 1)


@index_app.command("build")
def index_build(
    limit: Annotated[
        int | None,
        typer.Option(min=1, help="Download and index the first N PubChem SDF chunks."),
    ] = None,
    all_chunks: Annotated[
        bool, typer.Option("--all", help="Download and index all of PubChem.")
    ] = False,
    sdf: Annotated[
        list[Path] | None,
        typer.Option(
            "--sdf",
            exists=True,
            dir_okay=False,
            help="Index local SDF/SDF.gz files instead of downloading (repeatable).",
        ),
    ] = None,
    workers: Annotated[
        int, typer.Option(min=1, help="Parallel fingerprinting processes.")
    ] = _default_workers(),
    radius: Annotated[int, typer.Option(min=1, max=6, help="Morgan radius.")] = 2,
    fp_size: Annotated[
        int, typer.Option(min=64, max=8192, help="Fingerprint length in bits.")
    ] = 1024,
    keep_sdf: Annotated[
        bool, typer.Option(help="Keep downloaded SDF files after indexing them.")
    ] = False,
    home: HomeOpt = None,
) -> None:
    """Build (or extend) the local index, then merge it so it is searchable.

    Re-running is safe: finished chunks are skipped, so an interrupted build
    resumes and a larger --limit adds chunks to an existing index.
    """
    from .index.builder import SdfSource, build_chunks, finalize
    from .index.layout import IndexPaths, chunk_name

    chosen = sum([limit is not None, all_chunks, bool(sdf)])
    if chosen != 1:
        typer.echo(
            "error: pass exactly one of --limit N, --all or --sdf FILE", err=True
        )
        raise typer.Exit(code=2)

    paths = IndexPaths(resolve_home(home))
    params = {"radius": radius, "fpSize": fp_size}

    def say(msg: str) -> None:
        typer.echo(msg, err=True)

    try:
        if sdf:
            sources = [SdfSource(name=chunk_name(p), local_path=p) for p in sdf]
            report = build_chunks(paths, sources, workers, params, say)
        else:
            from .data.downloads import SdfMirror

            with SdfMirror() as mirror:
                names = mirror.list_chunks()
                if limit is not None:
                    names = names[:limit]
                say(
                    f"{len(names)} PubChem SDF chunk(s) selected; data directory {paths.home}"
                )

                def fetcher(name: str):  # type: ignore[no-untyped-def]
                    def fetch() -> Path:
                        say(f"downloading {name}")
                        return mirror.download(name, paths.downloads)

                    return fetch

                sources = [
                    SdfSource(
                        name=chunk_name(n), fetch=fetcher(n), delete_after=not keep_sdf
                    )
                    for n in names
                ]
                report = build_chunks(paths, sources, workers, params, say)
        manifest = finalize(paths, say)
    except PubChemFinderError as exc:
        _fail(exc)
        return
    typer.echo(
        f"Index ready: {manifest.n_compounds:,} compounds from "
        f"{len(manifest.chunks)} chunk(s), {len(report.built)} built this run, in {paths.home}"
    )


@index_app.command("status")
def index_status_cmd(home: HomeOpt = None, as_json: JsonFlag = False) -> None:
    """Show whether the local index is built and what it contains."""
    _emit(service.index_status(home), as_json)


@app.command()
def mcp() -> None:
    """Run the MCP server over stdio (for Claude Desktop, Claude Code, etc.)."""
    from .mcp_server import main

    main()


def run() -> None:  # pragma: no cover - console script entry
    sys.exit(app())


if __name__ == "__main__":  # pragma: no cover
    run()
