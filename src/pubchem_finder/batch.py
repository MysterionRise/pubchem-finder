"""Batch mode: run one operation over a file of structures and write a CSV.

Input formats:

- ``.csv`` / ``.tsv``: a header row with a SMILES column (``smiles``,
  ``canonical_smiles``, ``isomeric_smiles`` or ``smi``) and optionally an id
  column (``id``, ``name``, ``cid``, ``compound_id`` or ``title``).
- anything else (``.smi``, ``.txt``, ``-`` for stdin): one structure per line,
  ``SMILES [id ...]``; blank lines and ``#`` comments are skipped.

Rows that fail (invalid SMILES, no index, ...) are written with an ``error``
column instead of stopping the run.
"""

import csv
import sys
from collections.abc import Callable, Iterable, Iterator, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import IO, Any

from .core.exceptions import InvalidInput, PubChemFinderError

SMILES_COLUMNS = ("smiles", "canonical_smiles", "isomeric_smiles", "smi")
ID_COLUMNS = ("id", "name", "cid", "compound_id", "title")


@dataclass(frozen=True)
class InputRow:
    line: int
    id: str
    smiles: str


@dataclass
class BatchReport:
    rows_in: int = 0
    rows_out: int = 0
    errors: int = 0


def _open_input(path: Path) -> IO[str]:
    if str(path) == "-":
        return sys.stdin
    return open(path, newline="", encoding="utf-8")


def read_inputs(path: Path) -> Iterator[InputRow]:
    """Yield the structures in ``path`` (see module docstring for formats)."""
    suffix = path.suffix.lower()
    fh = _open_input(path)
    try:
        if suffix in (".csv", ".tsv"):
            yield from _read_table(fh, "\t" if suffix == ".tsv" else ",")
        else:
            yield from _read_smi(fh)
    finally:
        if fh is not sys.stdin:
            fh.close()


def _pick(headers: Sequence[str], names: Sequence[str]) -> str | None:
    lowered = {h.strip().lower(): h for h in headers}
    return next((lowered[n] for n in names if n in lowered), None)


def _read_table(fh: IO[str], delimiter: str) -> Iterator[InputRow]:
    reader = csv.DictReader(fh, delimiter=delimiter)
    headers = reader.fieldnames or []
    smiles_col = _pick(headers, SMILES_COLUMNS)
    if smiles_col is None:
        raise InvalidInput(
            f"No SMILES column found; expected one of {', '.join(SMILES_COLUMNS)}, "
            f"got {', '.join(headers) or 'no header'}"
        )
    id_col = _pick(headers, ID_COLUMNS)
    for n, record in enumerate(reader, start=2):
        smiles = (record.get(smiles_col) or "").strip()
        ident = (record.get(id_col) or "").strip() if id_col else ""
        yield InputRow(line=n, id=ident or str(n), smiles=smiles)


def _read_smi(fh: IO[str]) -> Iterator[InputRow]:
    for n, raw in enumerate(fh, start=1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split(maxsplit=1)
        ident = parts[1].strip() if len(parts) > 1 else str(n)
        yield InputRow(line=n, id=ident, smiles=parts[0])


def _cell(value: Any) -> Any:
    if value is None:
        return ""
    if isinstance(value, float):
        return round(value, 4)
    if isinstance(value, list | tuple):
        return "; ".join(str(v) for v in value)
    return value


def run(
    rows: Iterable[InputRow],
    out: IO[str],
    columns: Sequence[str],
    fn: Callable[[InputRow], list[dict[str, Any]]],
) -> BatchReport:
    """Apply ``fn`` to every row and stream the results to ``out`` as CSV.

    ``fn`` returns one or more output rows per input row. A
    ``PubChemFinderError`` becomes a single row with the ``error`` column set.
    """
    fields = list(dict.fromkeys(["id", "smiles", *columns, "error"]))
    writer = csv.DictWriter(out, fieldnames=fields, extrasaction="ignore")
    writer.writeheader()
    report = BatchReport()
    for row in rows:
        report.rows_in += 1
        base = {"id": row.id, "smiles": row.smiles}
        try:
            results = fn(row) if row.smiles else []
            if not row.smiles:
                raise InvalidInput(f"Empty SMILES on line {row.line}")
        except PubChemFinderError as exc:
            report.errors += 1
            results = [{"error": str(exc)}]
        for result in results or [{}]:
            writer.writerow({k: _cell(v) for k, v in {**base, **result}.items()})
            report.rows_out += 1
    return report
