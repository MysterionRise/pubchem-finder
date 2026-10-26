"""Command-line interface: ``pubchem-finder``."""

import os
import sys
from collections.abc import Callable
from enum import StrEnum
from pathlib import Path
from typing import Annotated, Any, cast

import typer
from pydantic import BaseModel

from . import __version__, batch, service
from .config import HOME_ENV, resolve_home
from .core.exceptions import InvalidInput, PubChemFinderError
from .data.pug_rest import ANNOTATION_HEADINGS, AnnotationHeading
from .models import RemoteSearchResults, SubstructureResults

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


class Transport(StrEnum):
    stdio = "stdio"
    http = "http"


JsonFlag = Annotated[bool, typer.Option("--json", help="Emit JSON instead of text.")]
SmilesArg = Annotated[
    str | None, typer.Argument(help="SMILES string (omit when using --input).")
]
InputOpt = Annotated[
    Path | None,
    typer.Option(
        "--input",
        "-i",
        help="Batch mode: .smi, .csv or .tsv file of structures ('-' for stdin).",
        show_default=False,
    ),
]
OutputOpt = Annotated[
    Path | None,
    typer.Option(
        "--output",
        "-o",
        help="Batch mode: CSV file to write (default: stdout).",
        show_default=False,
    ),
]


class HighlightMode(StrEnum):
    none = "none"
    smarts = "smarts"
    diff = "diff"
    similarity_map = "similarity_map"


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
    help="Build and inspect the local index (similarity and substructure search).",
    no_args_is_help=True,
    rich_markup_mode=None,
)
app.add_typer(index_app, name="index")

pubchem_app = typer.Typer(
    help="Query PubChem live: structure searches, synonyms, annotations, bioactivity.",
    no_args_is_help=True,
    rich_markup_mode=None,
)
app.add_typer(pubchem_app, name="pubchem")

cache_app = typer.Typer(
    help="Inspect or clear the cache of PubChem responses.",
    no_args_is_help=True,
    rich_markup_mode=None,
)
app.add_typer(cache_app, name="cache")


class Identity(StrEnum):
    same_stereo_isotope = "same_stereo_isotope"
    same_stereo = "same_stereo"
    same_isotope = "same_isotope"
    same_connectivity = "same_connectivity"
    same_tautomer = "same_tautomer"
    nonconflicting_stereo = "nonconflicting_stereo"


Heading = StrEnum(  # type: ignore[misc]
    "Heading", {h.replace(" ", "_").lower(): h for h in ANNOTATION_HEADINGS}
)


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
    smiles: SmilesArg = None,
    input_file: InputOpt = None,
    output: OutputOpt = None,
    as_json: JsonFlag = False,
) -> None:
    """RDKit descriptors and Lipinski check for a structure (or a file of them)."""
    if _batch_requested(smiles, input_file):
        _run_batch(
            input_file,
            output,
            PROPS_COLUMNS,
            lambda row: [
                service.properties(row.smiles).model_dump(exclude={"input_smiles"})
            ],
        )
        return
    try:
        _emit(service.properties(smiles or ""), as_json)
    except PubChemFinderError as exc:
        _fail(exc)


@app.command()
def standardize(
    smiles: SmilesArg = None,
    input_file: InputOpt = None,
    output: OutputOpt = None,
    as_json: JsonFlag = False,
) -> None:
    """Standardize a structure (remove salts, neutralize, canonicalize)."""
    if _batch_requested(smiles, input_file):

        def one(row: batch.InputRow) -> list[dict[str, Any]]:
            result = service.standardize(row.smiles)
            return [{"standardized_smiles": result.smiles, "changed": result.changed}]

        _run_batch(input_file, output, ["standardized_smiles", "changed"], one)
        return
    try:
        result = service.standardize(smiles or "")
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
    smiles: SmilesArg = None,
    threshold: Annotated[
        float, typer.Option(min=0.01, max=1.0, help="Minimum similarity (0-1).")
    ] = 0.7,
    top: Annotated[int, typer.Option(min=1, max=1000, help="Max hits.")] = 20,
    metric: Annotated[
        Metric, typer.Option(help="Similarity metric.")
    ] = Metric.tanimoto,
    group_by_scaffold: Annotated[
        bool, typer.Option(help="Also group the hits by Murcko scaffold.")
    ] = False,
    home: HomeOpt = None,
    input_file: InputOpt = None,
    output: OutputOpt = None,
    as_json: JsonFlag = False,
) -> None:
    """Find the most similar compounds in the local index (one query or a file)."""
    if _batch_requested(smiles, input_file):

        def one(row: batch.InputRow) -> list[dict[str, Any]]:
            found = service.similarity_search(
                row.smiles,
                threshold=threshold,
                top_k=top,
                metric=metric.value,
                home=home,
            )
            return [
                {"rank": rank, **hit.model_dump()}
                for rank, hit in enumerate(found.hits, start=1)
            ]

        _run_batch(input_file, output, SEARCH_COLUMNS, one)
        return
    try:
        result = service.similarity_search(
            smiles or "",
            threshold=threshold,
            top_k=top,
            metric=metric.value,
            home=home,
            group_by_scaffold=group_by_scaffold,
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
    if result.scaffolds is not None:
        typer.echo(f"\n{len(result.scaffolds)} scaffold(s):")
        for group in result.scaffolds:
            cids = ", ".join(str(c) for c in group.cids)
            typer.echo(
                f"  best {group.best_similarity:.3f}  {group.scaffold or '(acyclic)'}  CIDs {cids}"
            )


@app.command()
def substructure(
    query: Annotated[
        str | None,
        typer.Argument(
            help="SMILES (or SMARTS with --smarts) fragment; omit with --input."
        ),
    ] = None,
    smarts: Annotated[
        bool, typer.Option("--smarts", help="Read the query as SMARTS.")
    ] = False,
    max_hits: Annotated[
        int, typer.Option("--max", min=1, max=1000, help="Max hits.")
    ] = 50,
    max_candidates: Annotated[
        int,
        typer.Option(min=1, help="Stop after verifying this many screened candidates."),
    ] = 100_000,
    home: HomeOpt = None,
    input_file: InputOpt = None,
    output: OutputOpt = None,
    as_json: JsonFlag = False,
) -> None:
    """Find indexed compounds that contain a fragment, lowest CID first.

    Needs an index built with --with-substructure.
    """
    query_type: service.QueryType = "smarts" if smarts else "smiles"

    def run(text: str) -> SubstructureResults:
        return service.substructure_search(
            text,
            query_type=query_type,
            max_hits=max_hits,
            max_candidates=max_candidates,
            home=home,
        )

    if _batch_requested(query, input_file):

        def one(row: batch.InputRow) -> list[dict[str, Any]]:
            return [
                {"rank": rank, **hit.model_dump()}
                for rank, hit in enumerate(run(row.smiles).hits, start=1)
            ]

        _run_batch(input_file, output, SUBSTRUCTURE_COLUMNS, one)
        return
    try:
        result = run(query or "")
    except PubChemFinderError as exc:
        _fail(exc)
        return
    if as_json:
        typer.echo(result.model_dump_json(indent=2))
        return
    more = (
        ""
        if result.exhaustive
        else " (more may exist; raise --max or --max-candidates)"
    )
    typer.echo(
        f"{len(result.hits)} hit(s){more}; {result.n_candidates:,} of "
        f"{result.n_indexed:,} compounds passed the screen, {result.n_checked:,} checked"
    )
    if result.warning:
        typer.echo(f"warning: {result.warning}", err=True)
    for hit in result.hits:
        typer.echo(f"CID {hit.cid:<10} {hit.smiles or ''}  {hit.iupac_name or ''}")


@app.command()
def triage(
    smiles: SmilesArg = None,
    input_file: InputOpt = None,
    output: OutputOpt = None,
    as_json: JsonFlag = False,
) -> None:
    """Med-chem triage: QED, synthetic accessibility, Lipinski/Veber, structural alerts."""
    if _batch_requested(smiles, input_file):

        def one(row: batch.InputRow) -> list[dict[str, Any]]:
            profile = service.medchem_profile(row.smiles)
            data = profile.model_dump(exclude={"input_smiles", "alerts"})
            data["n_alerts"] = len(profile.alerts)
            data["alerts"] = [f"{a.catalog}:{a.name}" for a in profile.alerts]
            return [data]

        _run_batch(input_file, output, TRIAGE_COLUMNS, one)
        return
    try:
        result = service.medchem_profile(smiles or "")
    except PubChemFinderError as exc:
        _fail(exc)
        return
    if as_json:
        typer.echo(result.model_dump_json(indent=2))
        return
    sa = "-" if result.sa_score is None else f"{result.sa_score:.2f}"
    typer.echo(f"{result.canonical_smiles}")
    typer.echo(
        f"QED {result.qed:.2f}  SA score {sa}  scaffold {result.murcko_scaffold or '(acyclic)'}"
    )
    lipinski = ", ".join(result.lipinski_violations) or "none"
    typer.echo(
        f"Lipinski violations: {lipinski}   Veber: {'ok' if result.veber_ok else 'fails'}"
    )
    if not result.alerts:
        typer.echo("No PAINS, Brenk or NIH alerts.")
    for alert in result.alerts:
        typer.echo(f"  alert {alert.catalog}: {alert.name}  (atoms {alert.atoms})")


@app.command()
def draw(
    smiles: Annotated[list[str], typer.Argument(help="One or more SMILES (up to 24).")],
    output: Annotated[Path, typer.Option("--output", "-o", help="PNG file to write.")],
    legend: Annotated[
        list[str] | None,
        typer.Option("--legend", "-l", help="Caption per molecule (repeat, in order)."),
    ] = None,
    highlight: Annotated[
        HighlightMode, typer.Option(help="What to highlight.")
    ] = HighlightMode.none,
    smarts: Annotated[
        str | None, typer.Option(help="SMARTS pattern for --highlight smarts.")
    ] = None,
    reference: Annotated[
        str | None,
        typer.Option(help="Reference SMILES for --highlight diff / similarity_map."),
    ] = None,
    size: Annotated[
        int, typer.Option(min=150, max=600, help="Panel size in pixels.")
    ] = 300,
    per_row: Annotated[int, typer.Option(min=1, max=8, help="Panels per row.")] = 4,
) -> None:
    """Draw molecules to a PNG grid, optionally highlighting a pattern or changes."""
    try:
        png = service.draw(
            smiles,
            legends=legend,
            highlight=highlight.value,
            smarts=smarts,
            reference=reference,
            size=size,
            per_row=per_row,
        )
    except PubChemFinderError as exc:
        _fail(exc)
        return
    output.write_bytes(png)
    typer.echo(f"Wrote {output} ({len(smiles)} molecule(s), {len(png) / 1024:.0f} KB)")


PROPS_COLUMNS = [
    "canonical_smiles",
    "molecular_weight",
    "logp",
    "num_h_donors",
    "num_h_acceptors",
    "tpsa",
    "num_rotatable_bonds",
    "num_aromatic_rings",
    "num_heavy_atoms",
    "lipinski_ok",
]
TRIAGE_COLUMNS = [
    "canonical_smiles",
    "qed",
    "sa_score",
    "molecular_weight",
    "logp",
    "num_h_donors",
    "num_h_acceptors",
    "tpsa",
    "num_rotatable_bonds",
    "lipinski_violations",
    "veber_ok",
    "pains_free",
    "n_alerts",
    "alerts",
    "murcko_scaffold",
]
SUBSTRUCTURE_COLUMNS = [
    "rank",
    "cid",
    "smiles",
    "iupac_name",
    "molecular_formula",
    "molecular_weight",
    "pubchem_url",
]
SEARCH_COLUMNS = [
    "rank",
    "cid",
    "similarity",
    "smiles",
    "iupac_name",
    "molecular_formula",
    "molecular_weight",
    "pubchem_url",
]


def _batch_requested(smiles: str | None, input_file: Path | None) -> bool:
    if (smiles is None) == (input_file is None):
        typer.echo("error: give either a SMILES argument or --input FILE", err=True)
        raise typer.Exit(code=2)
    return input_file is not None


def _run_batch(
    input_file: Path | None,
    output: Path | None,
    columns: list[str],
    fn: Callable[[batch.InputRow], list[dict[str, Any]]],
) -> None:
    assert input_file is not None
    try:
        rows = batch.read_inputs(input_file)
        if output is None or str(output) == "-":
            report = batch.run(rows, sys.stdout, columns, fn)
            where = "stdout"
        else:
            with output.open("w", newline="", encoding="utf-8") as fh:
                report = batch.run(rows, fh, columns, fn)
            where = str(output)
    except (InvalidInput, OSError) as exc:
        _fail(exc)
        return
    typer.echo(
        f"{report.rows_in} input(s), {report.rows_out} row(s) written to {where}, "
        f"{report.errors} error(s)",
        err=True,
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
    with_substructure: Annotated[
        bool | None,
        typer.Option(
            "--with-substructure/--without-substructure",
            help=(
                "Add (or remove) substructure search. Once added it is kept up to "
                "date by later builds. Alone, adds it to the chunks already built."
            ),
            show_default=False,
        ),
    ] = None,
    home: HomeOpt = None,
) -> None:
    """Build (or extend) the local index, then merge it so it is searchable.

    Re-running is safe: finished chunks are skipped, so an interrupted build
    resumes and a larger --limit adds chunks to an existing index.
    """
    from .index.builder import BuildReport, SdfSource, build_chunks, finalize
    from .index.layout import IndexPaths, chunk_name

    chosen = sum([limit is not None, all_chunks, bool(sdf)])
    if chosen > 1 or (chosen == 0 and with_substructure is None):
        typer.echo(
            "error: pass exactly one of --limit N, --all or --sdf FILE "
            "(or only --with-substructure / --without-substructure)",
            err=True,
        )
        raise typer.Exit(code=2)

    paths = IndexPaths(resolve_home(home))
    params = {"radius": radius, "fpSize": fp_size}

    def say(msg: str) -> None:
        typer.echo(msg, err=True)

    try:
        if chosen == 0:
            report = BuildReport()
        elif sdf:
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
        manifest = finalize(paths, say, substructure=with_substructure, workers=workers)
    except PubChemFinderError as exc:
        _fail(exc)
        return
    extra = " with substructure search" if manifest.substructure else ""
    typer.echo(
        f"Index ready: {manifest.n_compounds:,} compounds from "
        f"{len(manifest.chunks)} chunk(s){extra}, {len(report.built)} built this run, "
        f"in {paths.home}"
    )


@index_app.command("status")
def index_status_cmd(home: HomeOpt = None, as_json: JsonFlag = False) -> None:
    """Show whether the local index is built and what it contains."""
    _emit(service.index_status(home), as_json)


def _emit_hits(result: RemoteSearchResults, as_json: bool) -> None:
    if as_json:
        typer.echo(result.model_dump_json(indent=2))
        return
    more = " (more exist; raise --max)" if result.truncated else ""
    typer.echo(f"{len(result.hits)} PubChem {result.search} hit(s){more}")
    for hit in result.hits:
        sim = "  -  " if hit.ecfp4_similarity is None else f"{hit.ecfp4_similarity:.3f}"
        typer.echo(f"{sim}  CID {hit.cid:<10} {hit.smiles or ''}  {hit.title or ''}")


@pubchem_app.command("similarity")
def pubchem_similarity(
    smiles: Annotated[str, typer.Argument(help="Query SMILES.")],
    threshold: Annotated[
        int, typer.Option(min=0, max=100, help="PubChem similarity threshold, percent.")
    ] = 90,
    max_records: Annotated[
        int, typer.Option("--max", min=1, max=200, help="Max hits.")
    ] = 25,
    as_json: JsonFlag = False,
) -> None:
    """Similar compounds across all of PubChem (ranked by ECFP4 Tanimoto)."""
    try:
        result = service.pubchem_similarity_search(
            smiles, threshold=threshold, max_records=max_records
        )
    except PubChemFinderError as exc:
        _fail(exc)
        return
    _emit_hits(result, as_json)


@pubchem_app.command("substructure")
def pubchem_substructure(
    smiles: Annotated[str, typer.Argument(help="Substructure SMILES.")],
    max_records: Annotated[
        int, typer.Option("--max", min=1, max=200, help="Max hits.")
    ] = 25,
    as_json: JsonFlag = False,
) -> None:
    """Compounds across all of PubChem that contain the given structure."""
    try:
        result = service.pubchem_substructure_search(smiles, max_records=max_records)
    except PubChemFinderError as exc:
        _fail(exc)
        return
    _emit_hits(result, as_json)


@pubchem_app.command("exact")
def pubchem_exact(
    smiles: Annotated[str, typer.Argument(help="Query SMILES.")],
    identity: Annotated[
        Identity, typer.Option(help="What counts as identical.")
    ] = Identity.same_stereo_isotope,
    as_json: JsonFlag = False,
) -> None:
    """The PubChem record(s) for exactly this structure."""
    try:
        result = service.pubchem_exact_match(smiles, identity_type=identity.value)
    except PubChemFinderError as exc:
        _fail(exc)
        return
    _emit_hits(result, as_json)


@pubchem_app.command("synonyms")
def pubchem_synonyms(
    cid: Annotated[int, typer.Argument(min=1, help="PubChem CID.")],
    limit: Annotated[int, typer.Option(min=1, max=500, help="Max names.")] = 25,
    as_json: JsonFlag = False,
) -> None:
    """Names and CAS numbers of a compound."""
    try:
        result = service.pubchem_synonyms(cid, limit=limit)
    except PubChemFinderError as exc:
        _fail(exc)
        return
    if as_json:
        typer.echo(result.model_dump_json(indent=2))
        return
    typer.echo(
        f"CAS: {', '.join(result.cas_numbers) or '-'}  ({result.total} names in PubChem)"
    )
    for name in result.synonyms:
        typer.echo(f"  {name}")


@pubchem_app.command("xrefs")
def pubchem_xrefs(
    cid: Annotated[int, typer.Argument(min=1, help="PubChem CID.")],
    as_json: JsonFlag = False,
) -> None:
    """Registry IDs, PubMed articles and patents linked to a compound."""
    try:
        result = service.pubchem_xrefs(cid)
    except PubChemFinderError as exc:
        _fail(exc)
        return
    if as_json:
        typer.echo(result.model_dump_json(indent=2))
        return
    for kind, count in result.counts.items():
        sample = ", ".join(result.samples[kind][:5])
        typer.echo(f"{kind:<11} {count:>7,}  {sample}")


@pubchem_app.command("describe")
def pubchem_describe(
    cid: Annotated[int, typer.Argument(min=1, help="PubChem CID.")],
    as_json: JsonFlag = False,
) -> None:
    """Curated descriptions of a compound, with sources."""
    try:
        result = service.pubchem_description(cid)
    except PubChemFinderError as exc:
        _fail(exc)
        return
    if as_json:
        typer.echo(result.model_dump_json(indent=2))
        return
    typer.echo(result.title or f"CID {cid}")
    for item in result.descriptions:
        typer.echo(f"\n{item.text}\n  -- {item.source or 'unknown source'}")


@pubchem_app.command("annotations")
def pubchem_annotations(
    cid: Annotated[int, typer.Argument(min=1, help="PubChem CID.")],
    heading: Annotated[Heading, typer.Option(help="Page section.")],
    max_chars: Annotated[
        int, typer.Option(min=500, max=30000, help="Text budget.")
    ] = 6000,
    as_json: JsonFlag = False,
) -> None:
    """One curated section of a compound's PubChem page (e.g. GHS safety)."""
    try:
        result = service.pubchem_annotations(
            cid, cast(AnnotationHeading, heading.value), max_chars=max_chars
        )
    except PubChemFinderError as exc:
        _fail(exc)
        return
    if as_json:
        typer.echo(result.model_dump_json(indent=2))
        return
    typer.echo(f"{result.title or cid}: {result.heading}")
    for a in result.annotations:
        label = f"{a.name}: " if a.name else ""
        typer.echo(f"- [{a.section}] {label}{a.text}  ({a.source or 'unknown source'})")
    if result.truncated:
        typer.echo(
            f"... {result.n_total - len(result.annotations)} more; raise --max-chars"
        )


@pubchem_app.command("bioactivity")
def pubchem_bioactivity(
    cid: Annotated[int, typer.Argument(min=1, help="PubChem CID.")],
    as_json: JsonFlag = False,
) -> None:
    """Summary of PubChem BioAssay results for a compound."""
    try:
        result = service.pubchem_bioactivity(cid)
    except PubChemFinderError as exc:
        _fail(exc)
        return
    if as_json:
        typer.echo(result.model_dump_json(indent=2))
        return
    if not result.tested:
        typer.echo(f"CID {cid} has no BioAssay results in PubChem.")
        return
    outcomes = ", ".join(f"{k} {v:,}" for k, v in result.outcomes.items())
    typer.echo(f"{result.n_assays:,} assays: {outcomes}")
    for t in result.top_targets:
        typer.echo(
            f"  target {t.target_accession or '-'} (gene {t.target_gene_id or '-'}):"
            f" active in {t.active_assays} assay(s)"
        )
    for p in result.most_potent:
        typer.echo(
            f"  {p.value_um:g} µM  AID {p.aid}  {p.activity_name or ''}  {p.assay_name or ''}"
        )
    typer.echo(result.pubchem_url)


@cache_app.command("info")
def cache_info(home: HomeOpt = None) -> None:
    """Show where the PubChem response cache is and how big it is."""
    stats = service.cache_stats(home)
    if stats is None:
        typer.echo("No PubChem response cache yet.")
        return
    typer.echo(
        f"{stats.path}: {stats.entries:,} responses, {stats.size_bytes / 1e6:.1f} MB"
    )


@cache_app.command("clear")
def cache_clear(home: HomeOpt = None) -> None:
    """Delete all cached PubChem responses."""
    typer.echo(f"Removed {service.cache_clear(home):,} cached responses.")


@app.command()
def mcp(
    transport: Annotated[
        Transport,
        typer.Option(
            "--transport",
            help="stdio: the client starts this process (default). "
            "http: serve MCP over streamable HTTP at http://HOST:PORT/mcp.",
        ),
    ] = Transport.stdio,
    host: Annotated[
        str,
        typer.Option(
            "--host",
            help="Interface to bind with --transport http; 127.0.0.1 keeps the server "
            "private to this machine.",
        ),
    ] = "127.0.0.1",
    port: Annotated[
        int,
        typer.Option(
            "--port", min=1, max=65535, help="TCP port with --transport http."
        ),
    ] = 8000,
    home: HomeOpt = None,
) -> None:
    """Run the MCP server over stdio (default, for Claude Desktop, Claude Code and
    other clients that start a local process) or over streamable HTTP."""
    if home is not None:
        os.environ[HOME_ENV] = str(resolve_home(home))
    from .mcp_server import main

    main(transport=transport.value, host=host, port=port)


def run() -> None:  # pragma: no cover - console script entry
    sys.exit(app())


if __name__ == "__main__":  # pragma: no cover
    run()
