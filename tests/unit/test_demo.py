"""The built-in demo index of FDA-approved drugs."""

import threading

import pytest
from typer.testing import CliRunner

from pubchem_finder import hints, service
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


# --- Searches fall back to the demo index ------------------------------------


@pytest.fixture
def demo_on(monkeypatch):
    monkeypatch.setenv("PUBCHEM_FINDER_DEMO", "1")


def test_similarity_search_without_an_index_searches_the_demo_set(demo_built, demo_on):
    found = service.similarity_search(ASPIRIN, threshold=0.4, home=demo_built)
    n = demo.provenance()["n_compounds"]
    assert found.search_scope == "demo"
    assert found.n_indexed == n
    assert found.hits[0].cid == 2244
    assert f"demo set of {n:,} FDA-approved drugs, not PubChem" in found.note
    assert "index build --limit 1 --with-substructure" in found.note
    assert "pubchem_similarity_search" in found.note
    # Through the search itself, not by formatting the note in the test: the
    # command is one that runs from any terminal, for the search that ran.
    assert f"`{hints.command('pubchem similarity SMILES')}`" in found.note
    assert "`pubchem-finder pubchem" not in found.note.replace(
        hints.command("pubchem similarity SMILES"), ""
    )


def test_substructure_search_without_an_index_searches_the_demo_set(
    demo_built, demo_on
):
    found = service.substructure_search("c1ccc2[nH]ccc2c1", home=demo_built)
    assert found.search_scope == "demo"
    assert found.hits
    assert f"`{hints.command('pubchem substructure SMILES')}`" in found.note
    assert "pubchem similarity SMILES" not in found.note


def test_status_says_searches_use_the_demo_set(demo_built, demo_on):
    status = service.index_status(demo_built)
    n = demo.provenance()["n_compounds"]
    assert (status.state, status.search_scope, status.demo_compounds) == (
        "not_built",
        "demo",
        n,
    )
    assert f"demo set of {n:,} FDA-approved drugs" in status.hint


def test_status_with_the_demo_off(tmp_path):
    status = service.index_status(tmp_path)
    assert (status.search_scope, status.demo_compounds) == ("none", 0)
    assert "demo set" not in status.hint


def test_a_built_index_takes_over_without_a_restart(tmp_path, sample_sdfs, demo_on):
    from pubchem_finder.index.builder import SdfSource, build_chunks, finalize
    from pubchem_finder.index.layout import chunk_name

    before = service.similarity_search(ASPIRIN, threshold=0.4, home=tmp_path)
    assert before.search_scope == "demo"
    paths = IndexPaths(tmp_path)
    build_chunks(paths, [SdfSource(chunk_name(p), local_path=p) for p in sample_sdfs])
    finalize(paths)
    after = service.similarity_search(ASPIRIN, threshold=0.4, home=tmp_path)
    assert (after.search_scope, after.note) == ("index", None)
    assert after.n_indexed == 11
    assert service.index_status(tmp_path).search_scope == "index"


def test_no_demo_fallback_for_substructure_on_an_index_without_it(built_home, demo_on):
    from pubchem_finder.core.exceptions import IndexNotBuilt

    with pytest.raises(IndexNotBuilt, match="no substructure search"):
        service.substructure_search("c1ccccc1", home=built_home)


def test_with_the_demo_off_searches_say_how_to_build_an_index(tmp_path):
    from pubchem_finder.core.exceptions import IndexNotBuilt

    with pytest.raises(IndexNotBuilt, match="index build --limit 1"):
        service.similarity_search(ASPIRIN, home=tmp_path)
    assert not (tmp_path / "demo").exists()


@pytest.mark.anyio
async def test_mcp_search_reports_the_demo_scope(demo_built, demo_on, monkeypatch):
    from mcp.client.client import Client

    from pubchem_finder.mcp_server import mcp

    monkeypatch.setenv("PUBCHEM_FINDER_HOME", str(demo_built))
    async with Client(mcp) as client:
        status = await client.call_tool("index_status", {})
        result = await client.call_tool(
            "similarity_search", {"smiles": ASPIRIN, "threshold": 0.4}
        )
    assert status.structured_content["search_scope"] == "demo"
    assert result.is_error is False
    assert result.structured_content["search_scope"] == "demo"
    assert result.structured_content["hits"][0]["cid"] == 2244
    assert "FDA-approved drugs" in result.structured_content["note"]


def test_cli_search_prints_the_demo_note(demo_built, demo_on):
    sim = runner.invoke(app, ["search", ASPIRIN, "--home", str(demo_built)])
    assert sim.exit_code == 0, sim.output
    assert "note: Searched the built-in demo set" in sim.output
    assert "CID 2244" in sim.output
    sub = runner.invoke(
        app, ["substructure", "c1ccc2[nH]ccc2c1", "--home", str(demo_built)]
    )
    assert sub.exit_code == 0, sub.output
    assert "note: Searched the built-in demo set" in sub.output
