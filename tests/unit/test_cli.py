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
