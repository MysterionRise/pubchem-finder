"""Tests for batch input parsing and CSV output."""

import csv
import io
from pathlib import Path

import pytest

from pubchem_finder import batch, service
from pubchem_finder.core.exceptions import InvalidInput

FIXTURES = Path(__file__).parents[1] / "fixtures" / "batch"


@pytest.mark.parametrize("name", ["utf8_bom_smiles_first.csv", "utf8_bom_id_first.csv"])
def test_read_csv_with_utf8_bom(name):
    rows = list(batch.read_inputs(FIXTURES / name))
    assert rows == [batch.InputRow(line=2, id="ethanol", smiles="CCO")]


@pytest.mark.parametrize(
    ("name", "text", "line"),
    [("in.tsv", "id\tsmiles\r\nethanol\tCCO\r\n", 2), ("in.smi", "CCO ethanol\n", 1)],
)
def test_read_other_formats_with_utf8_bom(tmp_path, name, text, line):
    path = tmp_path / name
    path.write_text(text, encoding="utf-8-sig")
    assert list(batch.read_inputs(path)) == [
        batch.InputRow(line=line, id="ethanol", smiles="CCO")
    ]


def test_read_non_utf8_input_has_actionable_error():
    path = FIXTURES / "windows_1252.csv"
    with pytest.raises(InvalidInput, match="UTF-8") as exc:
        list(batch.read_inputs(path))
    assert str(path) in str(exc.value)
    assert "save" in str(exc.value).lower()


def test_read_non_utf8_error_after_first_row(tmp_path):
    path = tmp_path / "late.csv"
    # Put the bad byte past the decoder's first buffer, after valid rows.
    path.write_bytes(b"id,smiles\n" + b"ethanol,CCO\n" * 1000 + b"caf\xe9,CCO\n")
    rows = batch.read_inputs(path)
    assert next(rows).smiles == "CCO"
    with pytest.raises(InvalidInput, match="UTF-8"):
        list(rows)


def write(tmp_path, name, text):
    path = tmp_path / name
    path.write_text(text)
    return path


def test_read_smi(tmp_path):
    path = write(
        tmp_path, "in.smi", "CCO ethanol\n\n# comment\nc1ccccc1\nCCN ethyl amine\n"
    )
    rows = list(batch.read_inputs(path))
    assert [(r.id, r.smiles) for r in rows] == [
        ("ethanol", "CCO"),
        ("4", "c1ccccc1"),
        ("ethyl amine", "CCN"),
    ]


def test_read_csv_and_tsv(tmp_path):
    csv_path = write(
        tmp_path,
        "in.csv",
        "Name,SMILES,other\naspirin,CC(=O)Oc1ccccc1C(=O)O,x\n,CCO,y\n",
    )
    rows = list(batch.read_inputs(csv_path))
    assert [(r.id, r.smiles) for r in rows] == [
        ("aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
        ("3", "CCO"),
    ]
    tsv_path = write(tmp_path, "in.tsv", "smi\nCCO\n")
    assert [r.smiles for r in batch.read_inputs(tsv_path)] == ["CCO"]


def test_csv_without_smiles_column(tmp_path):
    path = write(tmp_path, "in.csv", "name,structure\nx,CCO\n")
    with pytest.raises(InvalidInput, match="No SMILES column"):
        list(batch.read_inputs(path))


def test_read_stdin(monkeypatch):
    monkeypatch.setattr("sys.stdin", io.StringIO("CCO one\n"))
    assert [r.id for r in batch.read_inputs(batch.Path("-"))] == ["one"]


def test_run_writes_rows_and_errors():
    rows = [
        batch.InputRow(1, "a", "CCO"),
        batch.InputRow(2, "b", "nope("),
        batch.InputRow(3, "c", ""),
        batch.InputRow(4, "d", "c1ccccc1"),
    ]
    out = io.StringIO()

    def fn(row):
        props = service.properties(row.smiles)
        return [{"mw": props.molecular_weight, "tags": ["x", "y"], "none": None}]

    report = batch.run(rows, out, ["mw", "tags", "none"], fn)
    assert (report.rows_in, report.rows_out, report.errors) == (4, 4, 2)
    records = list(csv.DictReader(io.StringIO(out.getvalue())))
    assert [r["id"] for r in records] == ["a", "b", "c", "d"]
    assert records[0]["mw"] == "46.069"
    assert records[0]["tags"] == "x; y"
    assert records[0]["none"] == ""
    assert "Invalid SMILES" in records[1]["error"]
    assert "Empty SMILES on line 3" in records[2]["error"]


def test_run_multiple_and_zero_results_per_row():
    rows = [batch.InputRow(1, "a", "CCO"), batch.InputRow(2, "b", "CCN")]
    out = io.StringIO()
    report = batch.run(
        rows, out, ["rank"], lambda r: [{"rank": 1}, {"rank": 2}] if r.id == "a" else []
    )
    assert report.rows_out == 3  # two hits for a, one placeholder row for b
    records = list(csv.DictReader(io.StringIO(out.getvalue())))
    assert [(r["id"], r["rank"]) for r in records] == [
        ("a", "1"),
        ("a", "2"),
        ("b", ""),
    ]
