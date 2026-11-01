"""The built-in demo index of FDA-approved drugs."""

import threading

import pytest
from typer.testing import CliRunner

from pubchem_finder import service
from pubchem_finder.cli import app
from pubchem_finder.core.exceptions import IndexBuildError
from pubchem_finder.index import builder, demo
from pubchem_finder.index.layout import IndexPaths

ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"
runner = CliRunner()


@pytest.fixture(scope="module")
def demo_built(tmp_path_factory):
    """A data directory whose demo index is built (and nothing else)."""
    home = tmp_path_factory.mktemp("demo")
    demo.ensure_built(home)
    return home


def test_demo_index_holds_the_shipped_set_with_substructure_search(demo_built):
    manifest = demo.manifest(demo_built)
    assert manifest is not None
    assert manifest.n_compounds == demo.provenance()["n_compounds"]
    assert manifest.substructure is not None
    assert manifest.substructure.n_compounds == manifest.n_compounds
    assert demo.demo_home(demo_built).parent == demo_built / "demo"
    assert not IndexPaths(demo_built).manifest.exists()  # the main index is untouched


def test_demo_index_searches_like_any_index(demo_built):
    home = demo.demo_home(demo_built)
    found = service.similarity_search(ASPIRIN, threshold=0.4, home=home)
    assert found.hits[0].cid == 2244
    assert found.hits[0].iupac_name == "2-acetyloxybenzoic acid"
    indoles = service.substructure_search("c1ccc2[nH]ccc2c1", home=home, max_hits=5)
    assert len(indoles.hits) == 5


def test_built_demo_index_is_reused(demo_built, monkeypatch):
    def fail(*args, **kwargs):
        raise AssertionError("rebuilt")

    monkeypatch.setattr(demo, "_build", fail)
    assert demo.ensure_built(demo_built) == demo.demo_home(demo_built)


def test_key_changes_with_rdkit(monkeypatch):
    import rdkit

    before = demo.key()
    monkeypatch.setattr(rdkit, "__version__", "2000.01.1")
    assert demo.key() != before


def test_threads_build_it_once(tmp_path, monkeypatch):
    calls = []
    real = demo._build

    def counting(home, on_event):
        calls.append(home)
        real(home, on_event)

    monkeypatch.setattr(demo, "_build", counting)
    results = []
    threads = [
        threading.Thread(target=lambda: results.append(demo.ensure_built(tmp_path)))
        for _ in range(4)
    ]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    assert len(calls) == 1
    assert set(results) == {demo.demo_home(tmp_path)}


def test_losing_a_race_with_another_process_uses_its_index(tmp_path, monkeypatch):
    import os
    import shutil

    real_rename = os.rename

    def other_process_first(src, dst):
        shutil.copytree(src, dst)  # the other process's finished index
        real_rename(src, dst)  # now fails: the target exists

    monkeypatch.setattr(os, "rename", other_process_first)
    home = demo.ensure_built(tmp_path)
    assert IndexPaths(home).load_manifest() is not None
    assert [p.name for p in (tmp_path / "demo").iterdir()] == [home.name]


def test_failed_build_leaves_nothing_behind(tmp_path, monkeypatch):
    def broken(*args, **kwargs):
        raise OSError(28, "No space left on device")

    monkeypatch.setattr(demo, "write_chunk", broken)
    with pytest.raises(IndexBuildError, match="No space left"):
        demo.ensure_built(tmp_path)
    assert list((tmp_path / "demo").iterdir()) == []


def test_half_moved_index_is_replaced(tmp_path):
    target = demo.demo_home(tmp_path)
    (target / "chunks").mkdir(parents=True)  # no manifest: not a finished index
    assert demo.ensure_built(tmp_path) == target
    assert demo.manifest(tmp_path) is not None


def test_build_keeps_stdout_clean(tmp_path, capfd, monkeypatch):
    """Over stdio, stdout is the MCP transport."""
    real = builder.finalize

    def noisy(*args, **kwargs):
        print("FPSim2 says hello")
        return real(*args, **kwargs)

    monkeypatch.setattr(demo, "finalize", noisy)
    demo.ensure_built(tmp_path)
    out, err = capfd.readouterr()
    assert out == ""
    assert "FPSim2 says hello" in err


def test_tsv_reader_skips_unparsable_rows():
    text = (
        "cid\tsmiles\tinchikey\tiupac_name\tmolecular_formula\tmolecular_weight\n"
        "702\tCCO\t\t\t\t\n"
        "1\tC(C)(C)(C)(C)C\t\t\t\t1.0\n"
    )
    reader = demo.TsvReader(text)
    records = list(reader)
    assert [(r.cid, r.molecular_weight, r.iupac_name) for r in records] == [
        (702, None, None)
    ]
    assert (reader.n_ok, reader.n_failed) == (1, 1)


def test_cli_builds_the_demo_index(tmp_path):
    result = runner.invoke(app, ["index", "build", "--demo", "--home", str(tmp_path)])
    assert result.exit_code == 0, result.output
    n = demo.provenance()["n_compounds"]
    assert f"Demo index ready: {n:,} FDA-approved drugs" in result.output
    assert demo.manifest(tmp_path) is not None


@pytest.mark.parametrize("extra", [["--limit", "1"], ["--with-substructure"]])
def test_cli_demo_is_exclusive(tmp_path, extra):
    result = runner.invoke(
        app, ["index", "build", "--demo", *extra, "--home", str(tmp_path)]
    )
    assert result.exit_code == 2
    assert "--demo" in result.output


def test_cli_reports_a_failed_demo_build(tmp_path, monkeypatch):
    def broken(*args, **kwargs):
        raise OSError(13, "Permission denied")

    monkeypatch.setattr(demo, "write_chunk", broken)
    result = runner.invoke(app, ["index", "build", "--demo", "--home", str(tmp_path)])
    assert result.exit_code == 1
    assert "Permission denied" in result.output


def test_demo_switch(monkeypatch):
    from pubchem_finder.config import demo_enabled

    monkeypatch.delenv("PUBCHEM_FINDER_DEMO")
    assert demo_enabled()
    for off in ("0", "false", "OFF", " no "):
        monkeypatch.setenv("PUBCHEM_FINDER_DEMO", off)
        assert not demo_enabled()
