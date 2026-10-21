"""Tests for the typer CLI."""

import json

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
    assert "PubMedID" in xrefs.output and "CHEMBL25" in xrefs.output

    desc = runner.invoke(app, ["pubchem", "describe", "2244"])
    assert desc.output.startswith("Aspirin")
    assert "-- ChEBI" in desc.output

    ghs = runner.invoke(
        app, ["pubchem", "annotations", "2244", "--heading", "GHS Classification"]
    )
    assert ghs.exit_code == 0, ghs.output
    assert "Irritant" in ghs.output and "European Chemicals Agency" in ghs.output
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
    assert (
        "6 assays" in bio.output and "P23219" in bio.output and "0.5 µM" in bio.output
    )
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
