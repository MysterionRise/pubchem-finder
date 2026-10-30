"""Tests for the typer CLI."""

import csv
import json
import os

import pytest
from typer.testing import CliRunner

from pubchem_finder import __version__
from pubchem_finder.cli import app

runner = CliRunner()


def test_version():
    result = runner.invoke(app, ["--version"])
    assert result.exit_code == 0
    assert __version__ in result.output


def test_props_text_and_json(aspirin_smiles):
    text = runner.invoke(app, ["props", aspirin_smiles])
    assert text.exit_code == 0
    assert "molecular_weight" in text.output
    as_json = runner.invoke(app, ["props", aspirin_smiles, "--json"])
    assert as_json.exit_code == 0
    assert json.loads(as_json.output)["num_h_donors"] == 1


def test_props_invalid_smiles_exits_nonzero(invalid_smiles):
    result = runner.invoke(app, ["props", invalid_smiles])
    assert result.exit_code == 1
    assert "Invalid SMILES" in result.output


def test_standardize_prints_smiles():
    result = runner.invoke(app, ["standardize", "CC(=O)[O-].[Na+]"])
    assert result.exit_code == 0
    assert result.output.strip() == "CC(=O)O"


def test_compare():
    result = runner.invoke(app, ["compare", "CCO", "CCO", "--json"])
    assert result.exit_code == 0
    assert json.loads(result.output)["maccs_tanimoto"] == 1.0


def test_lookup_text_and_json(pubchem_api, monkeypatch):
    monkeypatch.setattr("pubchem_finder.data.pug_rest.time.sleep", lambda s: None)
    text = runner.invoke(app, ["lookup", "caffeine"])
    assert text.exit_code == 0
    assert "CID 2519  Caffeine" in text.output
    as_json = runner.invoke(app, ["lookup", "2244", "--by", "cid", "--json"])
    assert as_json.exit_code == 0
    assert json.loads(as_json.output)["compounds"][0]["title"] == "Aspirin"


def test_lookup_not_found(pubchem_api, monkeypatch):
    monkeypatch.setattr("pubchem_finder.data.pug_rest.time.sleep", lambda s: None)
    result = runner.invoke(app, ["lookup", "notacompound"])
    assert result.exit_code == 1
    assert "error:" in result.output


def test_index_build_requires_one_source(tmp_path):
    result = runner.invoke(app, ["index", "build", "--home", str(tmp_path)])
    assert result.exit_code == 2
    assert "exactly one" in result.output


def test_index_build_from_sdf_then_status_and_search(tmp_path, sample_sdfs):
    home = tmp_path / "home"
    args = ["index", "build", "--home", str(home), "--workers", "1"]
    for path in sample_sdfs[:2]:
        args += ["--sdf", str(path)]
    built = runner.invoke(app, args)
    assert built.exit_code == 0, built.output
    assert "Index ready: 11 compounds" in built.output

    status = runner.invoke(app, ["index", "status", "--home", str(home), "--json"])
    assert json.loads(status.output)["state"] == "ready"

    found = runner.invoke(
        app, ["search", "CC(=O)Oc1ccccc1C(=O)O", "--home", str(home), "--json"]
    )
    assert found.exit_code == 0
    assert json.loads(found.output)["hits"][0]["cid"] == 2244

    text = runner.invoke(
        app,
        ["search", "CC(=O)Oc1ccccc1C(=O)O", "--home", str(home), "--threshold", "0.3"],
    )
    assert "CID 2244" in text.output
    assert "2-acetyloxybenzoic acid" in text.output


def test_search_without_index(tmp_path):
    result = runner.invoke(app, ["search", "CCO", "--home", str(tmp_path)])
    assert result.exit_code == 1
    assert "index build" in result.output


def test_index_status_not_built(tmp_path):
    result = runner.invoke(app, ["index", "status", "--home", str(tmp_path)])
    assert result.exit_code == 0
    assert "not_built" in result.output


def test_index_build_explains_a_version_mismatch(tmp_path, sample_sdfs, monkeypatch):
    import rdkit

    home = ["--home", str(tmp_path), "--workers", "1"]
    with monkeypatch.context() as m:
        m.setattr(rdkit, "__version__", "2000.01.1")
        first = runner.invoke(
            app, ["index", "build", "--sdf", str(sample_sdfs[0]), *home]
        )
    assert first.exit_code == 0, first.output
    second = runner.invoke(app, ["index", "build", "--sdf", str(sample_sdfs[1]), *home])
    assert second.exit_code == 1
    assert "RDKit 2000.01.1" in second.output
    assert isinstance(second.exception, SystemExit)


def test_index_status_reports_a_damaged_manifest(tmp_path, built_home):
    import shutil

    home = tmp_path / "home"
    shutil.copytree(built_home, home)
    (home / "index" / "manifest.json").write_text("{")
    result = runner.invoke(app, ["index", "status", "--home", str(home)])
    assert result.exit_code == 1
    assert "error:" in result.output
    assert isinstance(result.exception, SystemExit)


def test_index_build_downloads_with_limit(tmp_path, sample_sdfs, monkeypatch):
    """--limit N lists the mirror, fetches the first N chunks and deletes them after."""
    import shutil

    from pubchem_finder.data import downloads

    fetched: list[str] = []

    class FakeMirror:
        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return None

        def list_chunks(self):
            return [p.name for p in sample_sdfs]

        def download(self, name, dest_dir):
            fetched.append(name)
            dest_dir.mkdir(parents=True, exist_ok=True)
            src = next(p for p in sample_sdfs if p.name == name)
            return shutil.copy(src, dest_dir / name)

    monkeypatch.setattr(downloads, "SdfMirror", FakeMirror)
    home = tmp_path / "home"
    result = runner.invoke(
        app, ["index", "build", "--limit", "1", "--home", str(home), "--workers", "1"]
    )
    assert result.exit_code == 0, result.output
    assert fetched == [sample_sdfs[0].name]
    assert "Index ready: 5 compounds" in result.output
    assert not list((home / "downloads").glob("*.sdf.gz"))  # cleaned up


def test_index_build_reports_errors(tmp_path, monkeypatch):
    from pubchem_finder.core.exceptions import DownloadError
    from pubchem_finder.data import downloads

    class BrokenMirror:
        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return None

        def list_chunks(self):
            raise DownloadError("mirror unreachable")

    monkeypatch.setattr(downloads, "SdfMirror", BrokenMirror)
    result = runner.invoke(app, ["index", "build", "--all", "--home", str(tmp_path)])
    assert result.exit_code == 1
    assert "mirror unreachable" in result.output


def test_pubchem_commands(pubchem_knowledge_api):
    aspirin = "CC(=O)Oc1ccccc1C(=O)O"
    sim = runner.invoke(app, ["pubchem", "similarity", aspirin, "--threshold", "90"])
    assert sim.exit_code == 0, sim.output
    assert "1.000  CID 2244" in sim.output
    assert "  -    CID 99999" in sim.output
    as_json = runner.invoke(
        app, ["pubchem", "similarity", aspirin, "--max", "4", "--json"]
    )
    assert json.loads(as_json.output)["truncated"] is True
    assert (
        "more exist"
        in runner.invoke(app, ["pubchem", "similarity", aspirin, "--max", "4"]).output
    )

    sub = runner.invoke(app, ["pubchem", "substructure", "c1ccccc1O"])
    assert "2 PubChem substructure hit(s)" in sub.output
    exact = runner.invoke(
        app, ["pubchem", "exact", aspirin, "--identity", "same_connectivity"]
    )
    assert "CID 2244" in exact.output

    syn = runner.invoke(app, ["pubchem", "synonyms", "2244", "--limit", "2"])
    assert "CAS: 50-78-2" in syn.output
    assert "ACETYLSALICYLIC ACID" in syn.output

    xrefs = runner.invoke(app, ["pubchem", "xrefs", "2244"])
    assert "PubMedID" in xrefs.output
    assert "CHEMBL25" in xrefs.output

    desc = runner.invoke(app, ["pubchem", "describe", "2244"])
    assert desc.output.startswith("Aspirin")
    assert "-- ChEBI" in desc.output

    ghs = runner.invoke(
        app, ["pubchem", "annotations", "2244", "--heading", "GHS Classification"]
    )
    assert ghs.exit_code == 0, ghs.output
    assert "Irritant" in ghs.output
    assert "European Chemicals Agency" in ghs.output
    short = runner.invoke(
        app,
        [
            "pubchem",
            "annotations",
            "2244",
            "--heading",
            "GHS Classification",
            "--max-chars",
            "500",
        ],
    )
    assert "more; raise --max-chars" in short.output

    bio = runner.invoke(app, ["pubchem", "bioactivity", "2244"])
    assert "6 assays" in bio.output
    assert "P23219" in bio.output
    assert "0.5 µM" in bio.output
    untested = runner.invoke(app, ["pubchem", "bioactivity", "5"])
    assert "no BioAssay results" in untested.output

    for args in (
        ["synonyms", "2244"],
        ["xrefs", "2244"],
        ["describe", "2244"],
        ["bioactivity", "2244"],
        ["annotations", "2244", "--heading", "GHS Classification"],
        ["substructure", "c1ccccc1O"],
        ["exact", aspirin],
    ):
        out = runner.invoke(app, ["pubchem", *args, "--json"])
        assert out.exit_code == 0, (args, out.output)
        json.loads(out.output)


def test_pubchem_command_errors(pubchem_knowledge_api):
    missing = runner.invoke(
        app, ["pubchem", "annotations", "2244", "--heading", "Toxicity"]
    )
    assert missing.exit_code == 1
    assert "no 'Toxicity' section" in missing.output
    bad = runner.invoke(app, ["pubchem", "similarity", "not smiles("])
    assert bad.exit_code == 1
    for args in (["synonyms", "5"], ["substructure", "not("], ["exact", "not("]):
        assert runner.invoke(app, ["pubchem", *args]).exit_code == 1
    for cmd in ("xrefs", "describe", "bioactivity"):
        assert (
            runner.invoke(app, ["pubchem", cmd, "0"]).exit_code == 2
        )  # CID must be >= 1


def test_cache_commands(tmp_path):
    from pubchem_finder.data.cache import ResponseCache

    empty = runner.invoke(app, ["cache", "info", "--home", str(tmp_path)])
    assert "No PubChem response cache" in empty.output
    ResponseCache(tmp_path / "cache" / "pubchem-http.sqlite").put("k", b"v")
    info = runner.invoke(app, ["cache", "info", "--home", str(tmp_path)])
    assert "1 responses" in info.output
    cleared = runner.invoke(app, ["cache", "clear", "--home", str(tmp_path)])
    assert "Removed 1 cached responses" in cleared.output


def test_triage_single_and_json():
    text = runner.invoke(app, ["triage", "O=C1C(=Cc2ccccc2)SC(=S)N1"])
    assert text.exit_code == 0, text.output
    assert "alert PAINS_A: ene_rhod_A(235)" in text.output
    clean = runner.invoke(app, ["triage", "CCO"])
    assert "No PAINS, Brenk or NIH alerts." in clean.output
    as_json = runner.invoke(app, ["triage", "CCO", "--json"])
    assert json.loads(as_json.output)["pains_free"] is True
    assert runner.invoke(app, ["triage", "nope("]).exit_code == 1


def test_draw_command(tmp_path):
    out = tmp_path / "grid.png"
    result = runner.invoke(
        app,
        [
            "draw",
            "CC(=O)Oc1ccccc1C(=O)O",
            "O=C(O)c1ccccc1O",
            "-l",
            "aspirin",
            "-l",
            "salicylic acid",
            "--highlight",
            "smarts",
            "--smarts",
            "[OX2H]",
            "-o",
            str(out),
        ],
    )
    assert result.exit_code == 0, result.output
    assert out.read_bytes().startswith(b"\x89PNG")
    bad = runner.invoke(app, ["draw", "CCO", "--highlight", "diff", "-o", str(out)])
    assert bad.exit_code == 1
    assert "reference" in bad.output


def test_batch_props_standardize_triage(tmp_path):
    src = tmp_path / "in.smi"
    src.write_text("CC(=O)[O-].[Na+] salt\nnope( broken\nCCO\n")
    for cmd in ("props", "standardize", "triage"):
        out = tmp_path / f"{cmd}.csv"
        result = runner.invoke(app, [cmd, "--input", str(src), "--output", str(out)])
        assert result.exit_code == 0, (cmd, result.output)
        assert "3 input(s), 3 row(s)" in result.output
        assert "1 error(s)" in result.output
        rows = list(csv.DictReader(out.open()))
        assert [r["id"] for r in rows] == ["salt", "broken", "3"]
        assert rows[1]["error"].startswith("Invalid SMILES")
    std = list(csv.DictReader((tmp_path / "standardize.csv").open()))
    assert std[0]["standardized_smiles"] == "CC(=O)O"
    assert std[0]["changed"] == "True"


def test_batch_to_stdout_and_argument_checks(tmp_path):
    src = tmp_path / "in.csv"
    src.write_text("id,smiles\nethanol,CCO\n")
    result = runner.invoke(app, ["props", "-i", str(src)])
    assert result.exit_code == 0
    assert "ethanol,CCO,CCO," in result.output
    both = runner.invoke(app, ["props", "CCO", "-i", str(src)])
    assert both.exit_code == 2
    neither = runner.invoke(app, ["props"])
    assert neither.exit_code == 2
    no_column = tmp_path / "bad.csv"
    no_column.write_text("name\nx\n")
    assert runner.invoke(app, ["props", "-i", str(no_column)]).exit_code == 1


def test_search_batch_and_scaffolds(tmp_path, built_home):
    src = tmp_path / "queries.smi"
    src.write_text("CC(=O)Oc1ccccc1C(=O)O aspirin\nCCO ethanol\n")
    out = tmp_path / "hits.csv"
    result = runner.invoke(
        app,
        [
            "search",
            "-i",
            str(src),
            "-o",
            str(out),
            "--home",
            str(built_home),
            "--threshold",
            "0.3",
        ],
    )
    assert result.exit_code == 0, result.output
    rows = list(csv.DictReader(out.open()))
    aspirin = [r for r in rows if r["id"] == "aspirin"]
    assert aspirin[0]["rank"] == "1"
    assert aspirin[0]["cid"] == "2244"
    assert any(r["id"] == "ethanol" for r in rows)
    grouped = runner.invoke(
        app,
        [
            "search",
            "CC(=O)Oc1ccccc1C(=O)O",
            "--home",
            str(built_home),
            "--threshold",
            "0.3",
            "--group-by-scaffold",
        ],
    )
    assert "scaffold(s):" in grouped.output
    assert "c1ccccc1" in grouped.output


def test_substructure_text_json_and_smarts(substructure_home):
    home = ["--home", str(substructure_home)]
    text = runner.invoke(app, ["substructure", "c1ccccc1C(=O)O", *home])
    assert text.exit_code == 0, text.output
    assert text.output.startswith("4 hit(s); ")
    assert "CID 2244" in text.output
    as_json = runner.invoke(
        app, ["substructure", "[OX2H]c", "--smarts", "--max", "1", "--json", *home]
    )
    data = json.loads(as_json.output)
    assert data["query_type"] == "smarts"
    assert [h["cid"] for h in data["hits"]] == [338]
    assert data["exhaustive"] is False
    capped = runner.invoke(
        app, ["substructure", "[#6]", "--smarts", "--max-candidates", "2", *home]
    )
    assert "more may exist" in capped.output
    bad = runner.invoke(app, ["substructure", "C1CC", *home])
    assert bad.exit_code == 1
    assert "Invalid SMILES" in bad.output


def test_substructure_batch(tmp_path, substructure_home):
    src = tmp_path / "fragments.smi"
    src.write_text("c1ccccc1C(=O)O benzoic\nC1CCCCC1 cyclohexane\nC1CC broken\n")
    out = tmp_path / "hits.csv"
    result = runner.invoke(
        app,
        [
            "substructure",
            "-i",
            str(src),
            "-o",
            str(out),
            "--home",
            str(substructure_home),
        ],
    )
    assert result.exit_code == 0, result.output
    assert "3 input(s)" in result.output
    assert "1 error(s)" in result.output
    rows = list(csv.DictReader(out.open()))
    benzoic = [r for r in rows if r["id"] == "benzoic"]
    assert [r["cid"] for r in benzoic] == ["243", "338", "2244", "4133"]
    assert benzoic[0]["rank"] == "1"
    assert [r["cid"] for r in rows if r["id"] == "cyclohexane"] == [""]
    assert next(r for r in rows if r["id"] == "broken")["error"]


def test_substructure_without_pattern_index(built_home):
    result = runner.invoke(app, ["substructure", "c1ccccc1", "--home", str(built_home)])
    assert result.exit_code == 1
    assert "--with-substructure" in result.output


def test_index_build_with_substructure_on_existing_chunks(tmp_path, sample_sdfs):
    home = ["--home", str(tmp_path / "home"), "--workers", "1"]
    built = runner.invoke(app, ["index", "build", "--sdf", str(sample_sdfs[0]), *home])
    assert built.exit_code == 0, built.output
    status = runner.invoke(app, ["index", "status", "--json", *home[:2]])
    assert json.loads(status.output)["substructure"] is False

    added = runner.invoke(app, ["index", "build", "--with-substructure", *home])
    assert added.exit_code == 0, added.output
    assert "with substructure search" in added.output
    assert "0 built this run" in added.output
    status = runner.invoke(app, ["index", "status", *home[:2]])
    assert "substructure         True" in status.output
    found = runner.invoke(app, ["substructure", "c1ccccc1", *home[:2]])
    assert "CID 243" in found.output

    removed = runner.invoke(app, ["index", "build", "--without-substructure", *home])
    assert removed.exit_code == 0, removed.output
    assert "substructure search removed" in removed.output


def test_index_build_with_substructure_needs_built_chunks(tmp_path):
    result = runner.invoke(
        app, ["index", "build", "--with-substructure", "--home", str(tmp_path)]
    )
    assert result.exit_code == 1
    assert "nothing to finalize" in result.output


def test_mcp_command_defaults_to_stdio(monkeypatch):
    calls: list[dict] = []
    monkeypatch.setattr("pubchem_finder.mcp_server.main", lambda **kw: calls.append(kw))
    result = runner.invoke(app, ["mcp"])
    assert result.exit_code == 0, result.output
    assert calls == [{"transport": "stdio", "host": "127.0.0.1", "port": 8000}]


def test_mcp_command_http_options(monkeypatch, tmp_path):
    calls: list[dict] = []
    monkeypatch.setattr("pubchem_finder.mcp_server.main", lambda **kw: calls.append(kw))
    monkeypatch.setenv("PUBCHEM_FINDER_HOME", str(tmp_path / "unused"))
    result = runner.invoke(
        app,
        [
            "mcp",
            "--transport",
            "http",
            "--host",
            "0.0.0.0",
            "--port",
            "9000",
            "--home",
            str(tmp_path),
        ],
    )
    assert result.exit_code == 0, result.output
    assert calls == [{"transport": "http", "host": "0.0.0.0", "port": 9000}]
    assert os.environ["PUBCHEM_FINDER_HOME"] == str(tmp_path.resolve())


@pytest.mark.parametrize(
    "args", [["--transport", "sse"], ["--port", "0"], ["--port", "70000"]]
)
def test_mcp_command_rejects_bad_options(monkeypatch, args):
    monkeypatch.setattr(
        "pubchem_finder.mcp_server.main",
        lambda **kw: pytest.fail("the server must not start"),
    )
    result = runner.invoke(app, ["mcp", *args])
    assert result.exit_code == 2
