"""Command-line interface: ``pubchem-finder``."""

import contextlib
import errno
import itertools
import os
import re
import secrets
import shutil
import stat
import sys
import tempfile
from collections.abc import Callable, Iterator
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path
from typing import IO, Annotated, Any, cast

import typer
from pydantic import BaseModel

from . import __version__, batch, service
from .config import HOME_ENV, resolve_home
from .core.exceptions import InvalidInput, PubChemFinderError
from .data.pug_rest import ANNOTATION_HEADINGS, AnnotationHeading, CompoundNotFound
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
        # It is written, not read. Asking for a readable path refused a
        # write-only file, and /dev/stdout on macOS, where a name for an open
        # descriptor can be read only if the descriptor can.
        readable=False,
    ),
]


class HighlightMode(StrEnum):
    none = "none"
    smarts = "smarts"
    diff = "diff"
    similarity_map = "similarity_map"


def _home_path(value: str) -> Path:
    if not value.strip():
        raise typer.BadParameter(
            "--home is empty; give a directory or leave the option out"
        )
    return Path(value)


HomeOpt = Annotated[
    Path | None,
    typer.Option(
        "--home",
        parser=_home_path,
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


def _emit(model: BaseModel, as_json: bool, lists: tuple[str, ...] = ()) -> None:
    """Print a model as JSON, or one "field value" line per field. A list is
    printed only when it is named in ``lists`` and has something in it."""
    if as_json:
        typer.echo(model.model_dump_json(indent=2))
        return
    for key, value in model.model_dump().items():
        if isinstance(value, list) and value and key in lists:
            value = ", ".join(str(item) for item in value)
        if value is None or isinstance(value, list):
            continue
        if isinstance(value, float):
            value = f"{value:.3f}"
        typer.echo(f"{key:<20} {value}")


def _warn_once(warning: str | None, seen: list[str]) -> None:
    """Batch mode: a warning goes to stderr the first time a row meets it; the
    output file has no place for it."""
    if warning and warning not in seen:
        seen.append(warning)
        typer.echo(f"warning: {warning}", err=True)


def _some(names: list[str]) -> str:
    more = f" and {len(names) - 3} more" if len(names) > 3 else ""
    return ", ".join(names[:3]) + more


def _fail(exc: Exception) -> None:
    # For a file error, show "<reason>: <path>" rather than str(exc): it drops the
    # "[Errno N]" noise and, unlike repr, keeps the path readable on Windows (an
    # OSError str doubles the backslashes there).
    if isinstance(exc, OSError) and exc.filename and exc.strerror:
        message = f"{exc.strerror}: {exc.filename}"
    else:
        message = str(exc)
    typer.echo(f"error: {message}", err=True)
    raise typer.Exit(code=1)


@app.command()
def lookup(
    query: Annotated[str, typer.Argument(help="Name, CID, SMILES or InChIKey.")],
    by: Annotated[By, typer.Option(help="Identifier type of QUERY.")] = By.name,
    limit: Annotated[int, typer.Option(min=1, help="Max records.")] = 5,
    home: HomeOpt = None,
    as_json: JsonFlag = False,
) -> None:
    """Look a compound up in PubChem."""
    try:
        result = service.lookup(query, by=by.value, limit=limit, home=home)
    except CompoundNotFound as exc:
        hint = (
            " QUERY was read as a name; for a CID, SMILES or InChIKey add "
            "--by cid, --by smiles or --by inchikey."
            if by is By.name
            else ""
        )
        _fail(PubChemFinderError(f"{exc}.{hint}" if hint else str(exc)))
        return
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
            return [
                {
                    "standardized_smiles": result.smiles,
                    "changed": result.changed,
                    "removed": ".".join(result.removed),
                    "note": result.note or "",
                }
            ]

        columns = ["standardized_smiles", "changed", "removed", "note"]
        _run_batch(input_file, output, columns, one)
        return
    try:
        result = service.standardize(smiles or "")
    except PubChemFinderError as exc:
        _fail(exc)
        return
    if as_json:
        typer.echo(result.model_dump_json(indent=2))
        return
    typer.echo(result.smiles)
    if result.removed:
        typer.echo(f"removed: {'.'.join(result.removed)}", err=True)
    if result.note:
        typer.echo(f"note: {result.note}", err=True)


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
    """Find the most similar compounds in the local index (one query or a file).

    With no index built, searches the built-in demo set of FDA-approved drugs.
    """
    if _batch_requested(smiles, input_file):
        warned: list[str] = []
        cut: list[str] = []

        def one(row: batch.InputRow) -> list[dict[str, Any]]:
            found = service.similarity_search(
                row.smiles,
                threshold=threshold,
                top_k=top,
                metric=metric.value,
                home=home,
            )
            _warn_once(found.warning, warned)
            if found.truncated:
                cut.append(row.id)
            return [
                {"rank": rank, **hit.model_dump()}
                for rank, hit in enumerate(found.hits, start=1)
            ]

        _run_batch(input_file, output, SEARCH_COLUMNS, one)
        if cut:
            # The CSV has no place for it, and a count of an input's rows would
            # be read as its number of hits at the threshold.
            more = (
                f"raise --top (up to {service.MAX_TOP_K}) or --threshold"
                if top < service.MAX_TOP_K
                else "raise --threshold"
            )
            typer.echo(
                f"note: {len(cut)} input(s) have more hits than --top {top} "
                f"({_some(cut)}): their rows are the best {top}, not all; {more} "
                "for complete lists",
                err=True,
            )
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
    except (PubChemFinderError, OSError) as exc:
        _fail(exc)
        return
    if as_json:
        typer.echo(result.model_dump_json(indent=2))
        return
    if result.note:
        typer.echo(f"note: {result.note}", err=True)
    if result.hint:
        typer.echo(f"hint: {result.hint}", err=True)
    if result.warning:
        typer.echo(f"warning: {result.warning}", err=True)
    # A list cut at --top is not the count at the threshold; say which it is.
    n = len(result.hits)
    count = f"top {n} of more than {n} hit(s)" if result.truncated else f"{n} hit(s)"
    typer.echo(
        f"{count} >= {result.threshold} {result.metric} "
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

    Needs an index built with --with-substructure; with no index built, searches
    the built-in demo set of FDA-approved drugs.
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
        warned: list[str] = []

        def one(row: batch.InputRow) -> list[dict[str, Any]]:
            found = run(row.smiles)
            _warn_once(found.warning, warned)
            return [
                {"rank": rank, **hit.model_dump()}
                for rank, hit in enumerate(found.hits, start=1)
            ]

        _run_batch(input_file, output, SUBSTRUCTURE_COLUMNS, one)
        return
    try:
        result = run(query or "")
    except (PubChemFinderError, OSError) as exc:
        _fail(exc)
        return
    if as_json:
        typer.echo(result.model_dump_json(indent=2))
        return
    if result.exhaustive:
        more = ""
    elif result.n_checked == result.n_candidates:
        # Every candidate was looked at, and some had no structure to check.
        more = " (more may exist; see the warning)"
    else:
        more = " (more may exist; raise --max or --max-candidates)"
    typer.echo(
        f"{len(result.hits)} hit(s){more}; {result.n_candidates:,} of "
        f"{result.n_indexed:,} compounds passed the screen, {result.n_checked:,} checked"
    )
    if result.note:
        typer.echo(f"note: {result.note}", err=True)
    if result.hint:
        typer.echo(f"hint: {result.hint}", err=True)
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
    try:
        output.write_bytes(png)
    except OSError as exc:
        _fail(exc)
        return
    typer.echo(f"Wrote {output} ({len(smiles)} molecule(s), {len(png) / 1024:.0f} KB)")


PROPS_COLUMNS = [
    "canonical_smiles",
    "molecular_formula",
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
    if input_file is None:  # callers check _batch_requested() first
        raise typer.BadParameter("an input file is required", param_hint="'--input'")
    if output is not None and str(output) != "-" and _same_file(input_file, output):
        raise typer.BadParameter(
            f"{output} is the input file; write the results to a new file",
            param_hint="'--output'",
        )
    try:
        reader = batch.read_inputs(input_file)
        # Read the first row now: that opens the input and checks its header before
        # anything is created for the output.
        first = next(reader, None)
        rows = itertools.chain([] if first is None else [first], reader)
        if output is None or str(output) == "-":
            report = batch.run(rows, sys.stdout, columns, fn)
            where = "stdout"
        else:
            with _replacing(output) as fh:
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
    if report.rows_in and report.errors == report.rows_in:
        raise typer.Exit(code=1)


# Paths that name an open descriptor rather than a file. `-i /dev/stdin` reads what
# `-i -` reads; `-o /dev/stdout`, or the /dev/fd/63 that `-o >(command)` becomes,
# writes to whatever that descriptor is: a pipe, or a file another program holds
# open, which must be written through the descriptor and not replaced.
_FD_PATH = re.compile(r"(?:/dev/fd|/proc/self/fd)/(\d+)")
_DESCRIPTOR_PATH = re.compile(
    r"/dev/(?:stdin|stdout|stderr|fd/\d+)|/proc/[^/]+/(?:task/[^/]+/)?fd/\d+"
)
_BINARY = getattr(os, "O_BINARY", 0)
_NOFOLLOW = getattr(os, "O_NOFOLLOW", 0)


def _input_descriptor(path: Path) -> int | None:
    """The open descriptor ``path`` stands for (0 for ``-``), or None for a file."""
    text = path.as_posix()
    if text in ("-", "/dev/stdin"):
        return 0
    named = _FD_PATH.fullmatch(text)
    return int(named.group(1)) if named else None


def _same_file(a: Path, b: Path) -> bool:
    """Whether the batch input ``a`` and the output ``b`` are one file.

    Standard input counts when it is redirected from a file (`< lib.smi`):
    reading a file while writing over it lost most of it. Piped input
    (`cat lib.smi | ...`) cannot be recognised; :func:`_replacing` keeps that
    case whole.
    """
    try:
        target = os.stat(b)
        descriptor = _input_descriptor(a)
        if descriptor is None:
            return os.path.samestat(os.stat(a), target)
        source = os.fstat(descriptor)
        # A terminal or a pipe is not a file that could be written over.
        return stat.S_ISREG(source.st_mode) and os.path.samestat(source, target)
    except (OSError, ValueError, OverflowError):
        # No such file or descriptor, or one that cannot be looked at. Nothing
        # says the two are one file, and opening them reports what is wrong.
        return False


def _names_descriptor(path: Path) -> bool:
    """Whether ``path`` is a name for an open descriptor, or a link to one
    (`ln -s /dev/stdout out.csv`, for a program that wants a file name)."""
    current = os.path.abspath(path)
    for _ in range(40):  # links to links, as far as the kernel follows them
        if _DESCRIPTOR_PATH.fullmatch(re.sub(r"^/+", "/", Path(current).as_posix())):
            return True
        try:
            link = os.readlink(current)
        except (OSError, ValueError):
            return False  # no link: an ordinary name, or nothing yet
        current = os.path.normpath(os.path.join(os.path.dirname(current), link))
    return False


def _replaceable(output: Path) -> tuple[Path, os.stat_result | None] | None:
    """The regular file a batch writing ``output`` replaces, with what is known of
    it now (None for a file that does not exist yet); or None when ``output`` has
    to be written the way :func:`open` writes it.

    That is the case for whatever is not a file of its own: a device, a FIFO, a
    name for an open descriptor. And for a link the kernel does not follow, be it
    dangling or refused (as one planted in a shared directory is): following it
    here would overrule the kernel.
    """
    if _names_descriptor(output):
        return None
    try:
        entry = os.lstat(output)
    except FileNotFoundError:
        return output, None  # or its directory is missing, which staging reports
    except OSError:
        return None
    if stat.S_ISREG(entry.st_mode):
        return output, entry
    if not stat.S_ISLNK(entry.st_mode):
        return None
    try:
        found = os.stat(output)  # the kernel follows the link, by its own rules
        resolved = Path(os.path.realpath(output))
        if stat.S_ISREG(found.st_mode) and os.path.samestat(found, os.stat(resolved)):
            return resolved, found  # the link stays, the file it names is replaced
    except OSError:
        pass
    return None


def _part_stem(name: str) -> str:
    """How the name of the partial file of an output called ``name`` begins: with
    all of it, unless the 14 characters that follow (``.<random>.part``) would
    make a name longer than a directory holds, 255 bytes."""
    stem = name
    while len(os.fsencode(stem)) > 255 - 14:
        stem = stem[:-1]
    return stem


@dataclass
class _Partial:
    """The file a batch's rows are staged in until the run ends."""

    handle: IO[str]
    path: Path
    beside: bool  # next to the output, where it can be renamed over it


def _stage_beside(target: Path, private: bool) -> _Partial:
    """Create the partial file for ``target`` next to it, under a name nothing has.

    No file is removed or written over to make room for it: a name that can be
    guessed could be the batch's own input, or a link someone put there. With
    ``private`` only this user can read it: an existing output's mode, group and
    access list say who may read its rows, and the partial file has none of them
    until it takes the output's place.
    """
    flags = os.O_RDWR | os.O_CREAT | os.O_EXCL | _BINARY
    while True:
        name = f"{_part_stem(target.name)}.{secrets.token_hex(4)}.part"
        part = target.with_name(name)
        try:
            descriptor = os.open(part, flags, 0o600 if private else 0o666)
        except FileExistsError:
            continue
        handle = os.fdopen(descriptor, "w", newline="", encoding="utf-8")
        return _Partial(handle, part, beside=True)


def _is(path: Path, known: os.stat_result) -> bool:
    """Whether ``path`` still names the file ``known`` describes."""
    try:
        now = os.lstat(path)
    except OSError:
        return False
    return (now.st_dev, now.st_ino) == (known.st_dev, known.st_ino)


def _reopen(path: Path, known: os.stat_result, flags: int) -> int:
    """Open ``path`` if it is still the file ``known`` describes.

    A batch can run for hours, in a directory where others may rename things.
    What lies under a name at the end is read or written only if it is the file
    this run made, or found at the start: never a link, never another file.
    """
    replaced = OSError(f"{path} was replaced by another file while the batch ran")
    try:
        descriptor = os.open(path, flags | _NOFOLLOW | _BINARY)
    except OSError as exc:
        if exc.errno == errno.ELOOP:  # a link now
            raise replaced from exc
        raise
    now = os.fstat(descriptor)
    if (now.st_dev, now.st_ino) != (known.st_dev, known.st_ino):
        os.close(descriptor)
        raise replaced
    return descriptor


def _attributes(file: int | Path) -> dict[str, bytes] | None:
    """A file's extended attributes (an access list, a security label) where the
    system lets Python read them, or None if they cannot be told."""
    if not hasattr(os, "listxattr"):
        return {}
    try:
        if isinstance(file, int):
            return {name: os.getxattr(file, name) for name in os.listxattr(file)}
        names = os.listxattr(file, follow_symlinks=False)
        return {n: os.getxattr(file, n, follow_symlinks=False) for n in names}
    except OSError as exc:
        return {} if exc.errno in (errno.ENOTSUP, errno.EOPNOTSUPP) else None


def _room(descriptor: int, folder: Path) -> int | None:
    """Free bytes on the filesystem an open file is on, or None where that cannot
    be told (a filesystem that reports no size at all reports no room either)."""
    try:
        if hasattr(os, "fstatvfs"):
            volume = os.fstatvfs(descriptor)
            return volume.f_bavail * volume.f_frsize if volume.f_blocks else None
        usage = shutil.disk_usage(folder)  # Windows: the folder's volume is the file's
        return usage.free if usage.total else None
    except OSError:
        return None


def _publish(
    partial: _Partial,
    mine: os.stat_result,
    target: Path,
    previous: os.stat_result | None,
) -> Path | None:
    """Make the finished partial file the content of ``target``.

    By renaming it into place where that changes nothing else: a new output, or
    one with content, a single name, this user's owner and group and no access
    list or attributes of its own. Any other output keeps what it has by being
    written in place: so is an empty one, which a caller may hold open to read
    the rows from (``tempfile.NamedTemporaryFile``), and one that cannot be
    replaced, such as a file mounted into a container or one in a directory
    that only lets its content change.

    ``mine`` identifies the partial file: see :func:`_reopen`. Returns the
    partial file if its rows were copied and it could then not be removed.
    """
    descriptor = partial.handle.fileno()
    # What the output is by now. Its owner may have changed its mode or group
    # while the batch ran, and the result is to have what the file has at the
    # end. If it was removed or moved away there is nothing to write into: the
    # rows become a new file of that name, with the mode the output had. What
    # was found at the start stands only where another file has taken the name.
    mode = None if previous is None else stat.S_IMODE(previous.st_mode)
    try:
        now = os.lstat(target)
    except FileNotFoundError:
        previous = None
    except OSError:
        pass
    else:
        same = previous is not None and (now.st_dev, now.st_ino) == (
            previous.st_dev,
            previous.st_ino,
        )
        if same:
            previous, mode = now, stat.S_IMODE(now.st_mode)
    rename = partial.beside and _is(partial.path, mine)
    if previous is not None:
        attributes = _attributes(descriptor)
        rename = (
            rename
            and previous.st_size > 0
            and previous.st_nlink == 1
            and (previous.st_uid, previous.st_gid) == (mine.st_uid, mine.st_gid)
            and attributes is not None
            and attributes == _attributes(target)
        )
    if rename and mode is not None:
        # Through the descriptor where the system can: the name may no longer
        # be this file's. (It cannot on Windows before Python 3.13, where an
        # open file keeps its name.)
        by_descriptor = os.chmod in os.supports_fd
        try:
            os.chmod(descriptor if by_descriptor else partial.path, mode)
        except OSError:
            rename = False
    partial.handle.close()  # Windows renames no open file
    if rename:
        try:
            os.replace(partial.path, target)
        except OSError:
            pass  # copied below, which asks less of the directory
        else:
            return None
    with contextlib.ExitStack() as opened:
        rows = opened.enter_context(
            os.fdopen(_reopen(partial.path, mine, os.O_RDONLY), "rb")
        )
        if previous is None:
            new = os.O_WRONLY | os.O_CREAT | os.O_EXCL | _NOFOLLOW | _BINARY
            created = os.open(target, new, 0o666 if mode is None else mode)
            out = opened.enter_context(os.fdopen(created, "wb"))
        else:
            out = opened.enter_context(
                os.fdopen(_reopen(target, previous, os.O_WRONLY), "wb")
            )
            # The old content goes before the new is in: while it is still there,
            # make sure the new will fit, on the filesystem the file is on (not
            # the folder's, for a file mounted into it).
            needed = mine.st_size - os.fstat(out.fileno()).st_size
            room = _room(out.fileno(), target.parent)
            if room is not None and needed > room:
                raise OSError(errno.ENOSPC, os.strerror(errno.ENOSPC))
            out.truncate(0)
        shutil.copyfileobj(rows, out)
    try:
        os.unlink(partial.path)
    except OSError:
        return partial.path  # a folder that lets files be made but not removed
    return None


def _stage(
    output: Path, target: Path, previous: os.stat_result | None
) -> _Partial | IO[str]:
    """Where the rows of a run go until it ends: a partial file beside the output,
    or in the system's temporary directory when nothing can be created beside an
    output that exists. Where neither can be had, the output itself, opened to be
    written directly: it is then not protected, and written as it always was.
    """
    try:
        if previous is None:
            try:
                return _stage_beside(target, private=False)
            except OSError:
                # If the output itself can be created, the trouble was the
                # partial file's (its longer name, say). If not, this is the
                # error to report.
                new = os.O_WRONLY | os.O_CREAT | os.O_EXCL | _BINARY
                descriptor = os.open(target, new, 0o666)
                return os.fdopen(descriptor, "w", newline="", encoding="utf-8")
        # Writable and not locked; its content is untouched. Not opened for
        # appending: an append-only file would pass, and can then be neither
        # replaced nor rewritten.
        os.close(os.open(target, os.O_WRONLY | _BINARY))
        with contextlib.suppress(OSError):
            return _stage_beside(target, private=True)
        try:
            descriptor, name = tempfile.mkstemp(
                prefix=f"{_part_stem(target.name)}.", suffix=".part"
            )
        except OSError:
            descriptor = _reopen(target, previous, os.O_WRONLY)
            os.ftruncate(descriptor, 0)
            return os.fdopen(descriptor, "w", newline="", encoding="utf-8")
        handle = os.fdopen(descriptor, "w", newline="", encoding="utf-8")
        return _Partial(handle, Path(name), beside=False)
    except OSError as exc:
        if exc.errno is None or exc.strerror is None:
            raise
        # Name the file that was asked for, not a partial file or a link's target.
        raise OSError(exc.errno, exc.strerror, str(output)) from exc


def _earlier_partial_files(target: Path, mine: Path) -> list[str]:
    stem = re.escape(_part_stem(target.name))
    partial = re.compile(stem + r"\.[0-9a-f]{8}\.part")
    try:
        with os.scandir(target.parent) as entries:
            names = [e.name for e in entries if partial.fullmatch(e.name)]
    except OSError:
        return []
    return sorted(name for name in names if name != mine.name)


def _from_interrupt(exc: BaseException | None) -> bool:
    """Whether Ctrl+C is behind ``exc``: itself, or wrapped by a library that was
    running Python code when it came (pydantic wraps one that arrives in a
    serializer in a ValueError)."""
    seen: set[int] = set()
    while exc is not None and id(exc) not in seen:
        if isinstance(exc, KeyboardInterrupt):
            return True
        seen.add(id(exc))
        exc = exc.__cause__ or exc.__context__
    return False


@contextlib.contextmanager
def _replacing(output: Path) -> Iterator[IO[str]]:
    """Open a text file whose content becomes ``output`` when the block succeeds.

    The rows are staged in a partial file that takes the output's place at the
    end. So a run that stops early leaves an existing output as it was, and input
    that is being piped from the output file is read to its end before the file
    changes. What cannot be written is reported before the batch runs, not after
    it: a missing directory, a read-only or locked output.

    After an error the partial file is removed. After Ctrl+C it is kept, and so
    it is when the run is killed: the rows finished so far may be hours of work.
    No run removes the partial file of another.

    A device, a pipe or an open descriptor (``-o /dev/stdout``) has nothing to
    keep and is written directly: see :func:`_replaceable`.
    """
    plan = _replaceable(output)
    if plan is None:
        with output.open("w", newline="", encoding="utf-8") as direct:
            yield direct
        return
    target, previous = plan
    staged = _stage(output, target, previous)
    if not isinstance(staged, _Partial):
        with staged as direct:
            yield direct
        return
    part = staged.path
    left = _earlier_partial_files(target, part) if staged.beside else []
    if left:
        typer.echo(
            f"note: {len(left)} partial file(s) of earlier runs lie beside {output} "
            f"({_some(left)}); this run leaves them as they are",
            err=True,
        )
    try:
        yield staged.handle
        staged.handle.flush()
        mine = os.fstat(staged.handle.fileno())
    except BaseException as exc:
        with contextlib.suppress(OSError):
            staged.handle.close()
        if _from_interrupt(exc):
            typer.echo(f"interrupted: the rows finished so far are in {part}", err=True)
            if isinstance(exc, KeyboardInterrupt):
                raise
            raise KeyboardInterrupt from exc
        with contextlib.suppress(OSError):  # the error that stopped the run counts
            part.unlink()
        if isinstance(exc, OSError) and exc.errno == errno.ENOSPC and not exc.filename:
            # A full disk names no file, and the rows may be staged on another
            # one than the output's.
            where = "" if staged.beside else " (set TMPDIR to a directory with room)"
            reason = f"{exc.strerror} where the rows are staged{where}"
            raise OSError(exc.errno, reason, str(part)) from exc
        raise
    try:
        copy = _publish(staged, mine, target, previous)
        if copy is not None:
            typer.echo(
                f"note: the rows are in {output}; the partial file they were "
                f"copied from could not be removed: {copy}",
                err=True,
            )
    except KeyboardInterrupt:
        if _is(part, mine):
            typer.echo(
                f"interrupted while {output} was being written: the complete rows "
                f"are in {part}",
                err=True,
            )
        raise
    except OSError as exc:
        kept = f"; the results are in {part}" if _is(part, mine) else ""
        raise OSError(
            f"could not write {output} ({exc.strerror or exc}){kept}"
        ) from exc
    finally:
        with contextlib.suppress(OSError):
            staged.handle.close()


def _default_workers() -> int:
    return max(1, (os.cpu_count() or 2) - 1)


def _download_progress(
    say: Callable[[str], None],
) -> Callable[[str, int, int | None], None]:
    """Report a download every tenth of the way (every 50 MB if the size is unknown)."""
    shown = [0]

    def report(name: str, done: int, total: int | None) -> None:
        step = done * 10 // total if total else done // 50_000_000
        if step > shown[0]:
            shown[0] = step
            size = f" of {total / 1e6:,.0f}" if total else ""
            say(f"  {done / 1e6:,.0f}{size} MB")

    return report


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
    demo: Annotated[
        bool,
        typer.Option(
            "--demo",
            help=(
                "Build the built-in demo index of FDA-approved drugs (a second, no "
                "download). Searches use it while no PubChem index is built."
            ),
        ),
    ] = False,
    workers: Annotated[
        int, typer.Option(min=1, help="Parallel fingerprinting processes.")
    ] = _default_workers(),
    radius: Annotated[
        int | None,
        typer.Option(
            min=1,
            max=6,
            help="Morgan radius (existing index setting, or 2).",
            show_default=False,
        ),
    ] = None,
    fp_size: Annotated[
        int | None,
        typer.Option(
            min=64,
            max=8192,
            help="Fingerprint length in bits (existing index setting, or 1024).",
            show_default=False,
        ),
    ] = None,
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
    from .index.builder import (
        DEFAULT_FP_PARAMS,
        BuildReport,
        SdfSource,
        build_chunks,
        finalize,
    )
    from .index.layout import IndexPaths, chunk_name

    chosen = sum([limit is not None, all_chunks, bool(sdf), demo])
    if (
        chosen > 1
        or (chosen == 0 and with_substructure is None)
        or (demo and with_substructure is not None)
    ):
        typer.echo(
            "error: pass exactly one of --limit N, --all, --sdf FILE or --demo "
            "(or only --with-substructure / --without-substructure)",
            err=True,
        )
        raise typer.Exit(code=2)

    if (demo or chosen == 0) and (radius is not None or fp_size is not None):
        typer.echo(
            "error: --radius / --fp-size applies only when chunks are built "
            "(--sdf, --limit, --all)",
            err=True,
        )
        raise typer.Exit(code=2)

    paths = IndexPaths(resolve_home(home))

    def say(msg: str) -> None:
        typer.echo(msg, err=True)

    if demo:
        _build_demo(paths.home, say)
        return

    try:
        completed = paths.completed_chunks()
        params = dict(completed[0].fp_params if completed else DEFAULT_FP_PARAMS)
        if radius is not None:
            params["radius"] = radius
        if fp_size is not None:
            params["fpSize"] = fp_size
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
                        return mirror.download(
                            name, paths.downloads, progress=_download_progress(say)
                        )

                    return fetch

                sources = [
                    SdfSource(
                        name=chunk_name(n), fetch=fetcher(n), delete_after=not keep_sdf
                    )
                    for n in names
                ]
                report = build_chunks(paths, sources, workers, params, say)
    except (PubChemFinderError, OSError) as exc:
        _merge_rebuilt_chunks(paths, say, workers)
        _fail(exc)
        return
    try:
        manifest = finalize(paths, say, substructure=with_substructure, workers=workers)
    except (PubChemFinderError, OSError) as exc:
        _fail(exc)
        return
    extra = " with substructure search" if manifest.substructure else ""
    typer.echo(
        f"Index ready: {manifest.n_compounds:,} compounds from "
        f"{len(manifest.chunks)} chunk(s){extra}, {len(report.built)} built this run, "
        f"in {paths.home}"
    )


def _merge_rebuilt_chunks(paths: Any, say: Callable[[str], None], workers: int) -> None:
    """After a build that failed part-way: merge the chunks it had replaced by then.

    Their new records are on disk while the index still holds the fingerprints of
    the files they replaced. Nothing is merged when the failed build replaced
    nothing: chunks that are only new wait for a build that succeeds.
    """
    from .index.builder import finalize, rebuilt_since_merge

    try:
        rebuilt = rebuilt_since_merge(paths)
        if not rebuilt:
            return
        say(
            f"{', '.join(rebuilt)}: indexed again before the build failed; merging "
            "so that the index matches the files"
        )
        finalize(paths, say, workers=workers)
    except (PubChemFinderError, OSError, ValueError) as exc:
        say(f"error: {exc}")


def _build_demo(home: Path, say: Callable[[str], None]) -> None:
    from .index import demo

    try:
        demo_home = demo.ensure_built(home, say)
    except (PubChemFinderError, OSError) as exc:
        _fail(exc)
        return
    info = demo.provenance()
    typer.echo(
        f"Demo index ready: {info['n_compounds']:,} {info['name']} with substructure "
        f"search, in {demo_home}"
    )


@index_app.command("status")
def index_status_cmd(home: HomeOpt = None, as_json: JsonFlag = False) -> None:
    """Show whether the local index is built and what it contains."""
    try:
        status = service.index_status(home)
    except (PubChemFinderError, OSError, ValueError) as exc:  # e.g. a damaged manifest
        _fail(exc)
        return
    _emit(status, as_json, lists=("unmerged_chunks",))


def _emit_hits(result: RemoteSearchResults, as_json: bool) -> None:
    if as_json:
        typer.echo(result.model_dump_json(indent=2))
        return
    more = " (more exist; raise --max)" if result.truncated else ""
    typer.echo(f"{len(result.hits)} PubChem {result.search} hit(s){more}")
    # Stderr keeps stdout's hit rows unchanged for scripts that parse them.
    if result.note:
        typer.echo(result.note, err=True)
    if result.hits:
        typer.echo(
            "score column: ECFP4 similarity to the query (computed locally)", err=True
        )
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
    home: HomeOpt = None,
    as_json: JsonFlag = False,
) -> None:
    """Similar compounds across all of PubChem (ranked by ECFP4 Tanimoto)."""
    try:
        result = service.pubchem_similarity_search(
            smiles, threshold=threshold, max_records=max_records, home=home
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
    home: HomeOpt = None,
    as_json: JsonFlag = False,
) -> None:
    """Compounds across all of PubChem that contain the given structure."""
    try:
        result = service.pubchem_substructure_search(
            smiles, max_records=max_records, home=home
        )
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
    max_records: Annotated[
        int, typer.Option("--max", min=1, max=200, help="Max hits.")
    ] = 25,
    home: HomeOpt = None,
    as_json: JsonFlag = False,
) -> None:
    """The PubChem record(s) for exactly this structure."""
    try:
        result = service.pubchem_exact_match(
            smiles, identity_type=identity.value, max_records=max_records, home=home
        )
    except PubChemFinderError as exc:
        _fail(exc)
        return
    _emit_hits(result, as_json)


@pubchem_app.command("synonyms")
def pubchem_synonyms(
    cid: Annotated[int, typer.Argument(min=1, help="PubChem CID.")],
    limit: Annotated[int, typer.Option(min=1, max=500, help="Max names.")] = 25,
    home: HomeOpt = None,
    as_json: JsonFlag = False,
) -> None:
    """Names and CAS numbers of a compound."""
    try:
        result = service.pubchem_synonyms(cid, limit=limit, home=home)
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
    home: HomeOpt = None,
    as_json: JsonFlag = False,
) -> None:
    """Registry IDs, PubMed articles and patents linked to a compound."""
    try:
        result = service.pubchem_xrefs(cid, home=home)
    except PubChemFinderError as exc:
        _fail(exc)
        return
    if as_json:
        typer.echo(result.model_dump_json(indent=2))
        return
    for kind, count in result.counts.items():
        sample = ", ".join(result.samples[kind][:5])
        typer.echo(f"{kind:<11} {count:>7,}  {sample}")
    for kind in result.unavailable:
        typer.echo(f"{kind:<11} {'?':>7}  (not available)")
    if result.note:
        typer.echo(result.note, err=True)


@pubchem_app.command("describe")
def pubchem_describe(
    cid: Annotated[int, typer.Argument(min=1, help="PubChem CID.")],
    home: HomeOpt = None,
    as_json: JsonFlag = False,
) -> None:
    """Curated descriptions of a compound, with sources."""
    try:
        result = service.pubchem_description(cid, home=home)
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
    home: HomeOpt = None,
    as_json: JsonFlag = False,
) -> None:
    """One curated section of a compound's PubChem page (e.g. GHS safety)."""
    try:
        result = service.pubchem_annotations(
            cid, cast(AnnotationHeading, heading.value), max_chars=max_chars, home=home
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
    home: HomeOpt = None,
    as_json: JsonFlag = False,
) -> None:
    """Summary of PubChem BioAssay results for a compound."""
    try:
        result = service.pubchem_bioactivity(cid, home=home)
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
