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


# Hardening audit D1/B5: chunks are named after their file, so vendor exports that
# share a name overwrote each other, and a changed file was skipped as built.


def _sdf_build(home, *files):
    args = ["index", "build", "--home", str(home), "--workers", "1"]
    for f in files:
        args += ["--sdf", str(f)]
    return runner.invoke(app, args)


def _cids(home):
    found = runner.invoke(
        app,
        ["search", "c1ccccc1", "--home", str(home), "--threshold", "0.01", "--json"],
    )
    return {hit["cid"] for hit in json.loads(found.output)["hits"]}


def test_build_says_when_ids_are_why_records_were_skipped(tmp_path):
    """A ChEMBL or vendor file (IDs like CHEMBL25) built 0 compounds and the only
    explanation was "No built chunks"."""
    from rdkit import Chem

    lib = tmp_path / "chembl.sdf"
    with Chem.SDWriter(str(lib)) as writer:
        for i, smiles in enumerate(["CCO", "CCN"], 1):
            mol = Chem.MolFromSmiles(smiles)
            mol.SetProp("_Name", f"CHEMBL{i}")
            writer.write(mol)
    result = _sdf_build(tmp_path / "home", lib)
    assert result.exit_code == 1
    assert "built chembl: 0 compounds (2 skipped)" in result.output
    assert "record IDs must be whole numbers from 1 to 4,294,967,295" in result.output
    assert "skipped 'CHEMBL1' and any others like it" in result.output


@pytest.fixture
def same_named(tmp_path, sample_sdfs):
    import shutil

    a = tmp_path / "vendor_a" / "library.sdf.gz"
    b = tmp_path / "vendor_b" / "library.sdf.gz"
    for path, src in ((a, sample_sdfs[0]), (b, sample_sdfs[1])):
        path.parent.mkdir()
        shutil.copyfile(src, path)
    return a, b


def test_index_build_refuses_two_files_with_one_name(tmp_path, same_named):
    result = _sdf_build(tmp_path / "home", *same_named)
    assert result.exit_code == 1
    assert "would both be chunk 'library'" in result.output
    assert not (tmp_path / "home" / "chunks").exists()


def test_index_build_refuses_a_file_that_would_replace_another(tmp_path, same_named):
    home = tmp_path / "home"
    assert _sdf_build(home, same_named[0]).exit_code == 0
    result = _sdf_build(home, same_named[1])
    assert result.exit_code == 1
    assert "already holds" in result.output
    assert 2244 in _cids(home)  # vendor A's compounds are still there


def test_index_build_rebuilds_a_changed_file(tmp_path, same_named, sample_sdfs):
    import shutil

    home = tmp_path / "home"
    library = same_named[0]
    assert _sdf_build(home, library).exit_code == 0
    again = _sdf_build(home, library)
    assert "already built, skipping" in again.output
    assert "0 built this run" in again.output

    shutil.copyfile(sample_sdfs[1], library)  # the vendor ships a new export
    result = _sdf_build(home, library)
    assert result.exit_code == 0, result.output
    assert "rebuilding library: " in result.output
    assert "changed since it was indexed" in result.output
    assert "Index ready: 6 compounds" in result.output
    cids = _cids(home)
    assert 1983 in cids  # from the new file
    assert 2244 not in cids  # gone with the old one


def _as_2_2_0_wrote_it(home, name="library"):
    stats = home / "chunks" / f"{name}.json"
    data = json.loads(stats.read_text())
    del data["source"], data["build_id"]
    stats.write_text(json.dumps(data))
    return stats


def test_index_build_adopts_the_file_of_a_chunk_from_an_earlier_version(
    tmp_path, same_named
):
    home = tmp_path / "home"
    library = same_named[0]
    assert _sdf_build(home, library).exit_code == 0
    stats = _as_2_2_0_wrote_it(home)
    result = _sdf_build(home, library)
    assert result.exit_code == 0, result.output
    assert "recorded" in result.output
    assert "0 built this run" in result.output  # not fingerprinted again
    assert json.loads(stats.read_text())["source"]["path"] == str(library.resolve())
    assert 2244 in _cids(home)


def test_index_build_rebuilds_an_earlier_versions_chunk_from_a_newer_file(
    tmp_path, same_named, sample_sdfs
):
    import shutil

    home = tmp_path / "home"
    library = same_named[0]
    assert _sdf_build(home, library).exit_code == 0
    stats = _as_2_2_0_wrote_it(home)
    shutil.copyfile(sample_sdfs[1], library)
    later = stats.stat().st_mtime_ns + 5_000_000_000
    os.utime(library, ns=(later, later))
    result = _sdf_build(home, library)
    assert result.exit_code == 0, result.output
    assert "is newer than its chunk" in result.output
    cids = _cids(home)
    assert 1983 in cids
    assert 2244 not in cids


# The 2.2.1 rework of --sdf rebuilds, reviewed after release: the old chunk was deleted
# before every check had run, and whether to merge was decided from what the process
# had just done, so a rebuild that failed part-way left an index that the next run
# called ready while it still searched the old file.


@pytest.mark.parametrize("fails_in", ["build_chunk", "_merge_h5"])
def test_index_build_after_a_failed_rebuild_serves_the_new_file(
    tmp_path, same_named, sample_sdfs, monkeypatch, fails_in
):
    import shutil

    from pubchem_finder.index import builder

    home = tmp_path / "home"
    library = same_named[0]
    assert _sdf_build(home, library).exit_code == 0
    shutil.copyfile(sample_sdfs[1], library)

    real = getattr(builder, fails_in)
    calls = []

    def fail_once(*args, **kwargs):
        calls.append(1)
        if len(calls) == 1:
            raise OSError("No space left on device")
        return real(*args, **kwargs)

    monkeypatch.setattr(builder, fails_in, fail_once)
    failed = _sdf_build(home, library)
    assert failed.exit_code == 1
    assert "No space left" in failed.output
    if fails_in == "build_chunk":  # nothing was replaced: the old index still answers
        assert 2244 in _cids(home)

    again = _sdf_build(home, library)
    assert again.exit_code == 0, again.output
    cids = _cids(home)
    assert 1983 in cids  # the new file
    assert 2244 not in cids  # not the old one


def test_refused_build_leaves_the_changed_files_chunk_alone(
    tmp_path, sample_sdfs, monkeypatch
):
    import shutil

    from pubchem_finder.index import builder

    home = tmp_path / "home"
    a = tmp_path / "a.sdf.gz"
    b = tmp_path / "b.sdf.gz"
    shutil.copyfile(sample_sdfs[0], a)
    shutil.copyfile(sample_sdfs[1], b)
    assert _sdf_build(home, a, b).exit_code == 0
    before = {p.name: p.stat().st_size for p in (home / "chunks").iterdir()}

    shutil.copyfile(sample_sdfs[0], b)  # b changes; a was built by "another RDKit"
    monkeypatch.setattr(builder, "runtime_versions", lambda: ("2099.01.1", "9.9.9"))
    refused = _sdf_build(home, b)
    assert refused.exit_code == 1
    assert "fingerprinted with" in refused.output
    assert {p.name: p.stat().st_size for p in (home / "chunks").iterdir()} == before
    monkeypatch.undo()
    assert {2244, 1983} <= _cids(home)  # both files' compounds, as before


def test_index_build_keeps_the_chunk_of_a_moved_file(tmp_path, same_named):
    import shutil

    home = tmp_path / "home"
    library = same_named[0]
    assert _sdf_build(home, library).exit_code == 0
    moved = tmp_path / "elsewhere" / library.name
    moved.parent.mkdir()
    shutil.move(library, moved)  # keeps size and modification time
    result = _sdf_build(home, moved)
    assert result.exit_code == 0, result.output
    assert "its file moved to" in result.output
    assert "0 built this run" in result.output
    stats = json.loads((home / "chunks" / "library.json").read_text())
    assert stats["source"]["path"] == str(moved.resolve())


def test_index_build_compares_chunk_names_without_case(tmp_path, sample_sdfs):
    import shutil

    upper = tmp_path / "vendor_a" / "Library.sdf.gz"
    lower = tmp_path / "vendor_b" / "library.sdf.gz"
    for path, src in ((upper, sample_sdfs[0]), (lower, sample_sdfs[1])):
        path.parent.mkdir()
        shutil.copyfile(src, path)
    both = _sdf_build(tmp_path / "home", upper, lower)
    assert both.exit_code == 1
    assert "without case" in both.output
    assert not (tmp_path / "home" / "chunks").exists()

    home = tmp_path / "home2"
    assert _sdf_build(home, upper).exit_code == 0
    later = _sdf_build(home, lower)
    assert later.exit_code == 1
    assert "already has chunk 'Library'" in later.output
    assert 2244 in _cids(home)


def test_index_build_knows_its_file_under_another_spelling(tmp_path, same_named):
    """On macOS and Windows a file can be typed in any case; a hard link stands in
    for that here, so the test means the same on every filesystem."""
    home = tmp_path / "home"
    library = same_named[0]
    assert _sdf_build(home, library).exit_code == 0
    alias = tmp_path / "alias" / "LIBRARY.sdf.gz"
    alias.parent.mkdir()
    os.link(library, alias)
    result = _sdf_build(home, alias)
    assert result.exit_code == 0, result.output
    assert "already built, skipping" in result.output
    assert sorted(p.name for p in (home / "chunks").glob("*.json")) == ["library.json"]


def test_index_build_says_when_a_file_is_given_twice(tmp_path, same_named):
    result = _sdf_build(tmp_path / "home", same_named[0], same_named[0])
    assert result.exit_code == 1
    assert "is given twice" in result.output


def test_index_build_rebuilds_a_changed_file_with_two_workers(
    tmp_path, same_named, sample_sdfs
):
    import shutil

    home = tmp_path / "home"
    library, other = same_named[0], tmp_path / "other.sdf.gz"
    shutil.copyfile(sample_sdfs[1], other)
    args = ["index", "build", "--home", str(home), "--workers", "2"]
    args += ["--sdf", str(library), "--sdf", str(other)]
    assert runner.invoke(app, args).exit_code == 0
    shutil.copyfile(sample_sdfs[1], library)
    result = runner.invoke(app, args)
    assert result.exit_code == 0, result.output
    assert "rebuilding library: " in result.output
    cids = _cids(home)
    assert 1983 in cids
    assert 2244 not in cids


def test_an_index_merged_by_2_2_1_with_sdf_chunks_is_merged_once_more(
    tmp_path, same_named, sample_sdfs
):
    """2.2.1 could rebuild a chunk and fail before merging it, leaving no trace in
    the manifest. Such an index is merged again the first time it is built on."""
    import shutil

    from pubchem_finder.index.builder import SdfSource, build_chunks
    from pubchem_finder.index.layout import IndexPaths

    home = tmp_path / "home"
    library = same_named[0]
    assert _sdf_build(home, library).exit_code == 0
    shutil.copyfile(sample_sdfs[1], library)
    build_chunks(IndexPaths(home), [SdfSource(name="library", local_path=library)])
    for path, field in (
        (home / "index" / "manifest.json", "chunk_builds"),
        (home / "chunks" / "library.json", "build_id"),
    ):  # as 2.2.1 left it: chunk rebuilt, nothing merged, no record of builds
        data = json.loads(path.read_text())
        del data[field]
        path.write_text(json.dumps(data))
    assert 2244 in _cids(home)  # stale: still the old file's fingerprints

    result = _sdf_build(home, library)
    assert result.exit_code == 0, result.output
    assert "merging" in result.output
    cids = _cids(home)
    assert 1983 in cids
    assert 2244 not in cids
    manifest = json.loads((home / "index" / "manifest.json").read_text())
    assert manifest["chunk_builds"] == {"library": None}
    assert "index up to date" in _sdf_build(home, library).output


def test_index_build_skips_ids_the_index_cannot_hold(tmp_path):
    from tests.conftest import write_sdf

    aspirin = "CC(=O)Oc1ccccc1C(=O)O"
    library = write_sdf(
        tmp_path / "registry.sdf.gz",
        [(2**32 + 2244, aspirin, "in-house number", "new"), (702, "CCO", "x", "new")],
    )
    home = tmp_path / "home"
    built = _sdf_build(home, library)
    assert built.exit_code == 0, built.output
    assert "1 compounds (1 skipped)" in built.output
    assert "record IDs must be whole numbers from 1 to 4,294,967,295" in built.output
    assert str(2**32 + 2244) in built.output
    found = runner.invoke(
        app, ["search", aspirin, "--home", str(home), "--threshold", "0.9", "--json"]
    )
    # 2.2.1 answered CID 2244 here, at similarity 1.0: the ID modulo 2**32.
    assert json.loads(found.output)["hits"] == []


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


def test_index_status_explains_a_newer_index_format(tmp_path, built_home):
    import json
    import shutil

    home = tmp_path / "home"
    shutil.copytree(built_home, home)
    manifest_file = home / "index" / "manifest.json"
    manifest = json.loads(manifest_file.read_text())
    manifest["format_version"] += 1
    manifest_file.write_text(json.dumps(manifest))
    result = runner.invoke(app, ["index", "status", "--home", str(home)])
    assert result.exit_code == 1
    assert "Upgrade pubchem-finder" in result.output
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

        def download(self, name, dest_dir, progress=None):
            fetched.append(name)
            dest_dir.mkdir(parents=True, exist_ok=True)
            src = next(p for p in sample_sdfs if p.name == name)
            for done in range(0, 323_000_001, 1_000_000):
                progress(name, done, 323_000_000)
            return shutil.copy(src, dest_dir / name)

    monkeypatch.setattr(downloads, "SdfMirror", FakeMirror)
    home = tmp_path / "home"
    result = runner.invoke(
        app, ["index", "build", "--limit", "1", "--home", str(home), "--workers", "1"]
    )
    assert result.exit_code == 0, result.output
    assert fetched == [sample_sdfs[0].name]
    progress = [line for line in result.output.splitlines() if line.endswith(" MB")]
    assert progress[0] == "  33 of 323 MB"
    assert progress[-1] == "  323 of 323 MB"
    assert len(progress) == 10
    assert "Index ready: 5 compounds" in result.output
    assert not list((home / "downloads").glob("*.sdf.gz"))  # cleaned up


def test_download_progress_without_a_size():
    from pubchem_finder.cli import _download_progress

    lines: list[str] = []
    report = _download_progress(lines.append)
    for done in range(0, 120_000_001, 10_000_000):
        report("Compound_x.sdf.gz", done, None)
    assert lines == ["  50 MB", "  100 MB"]


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


def test_pubchem_xrefs_says_which_types_are_missing(pubchem_knowledge_api):
    import httpx

    pubchem_knowledge_api["xrefs_patents"].mock(
        return_value=httpx.Response(504, text="Gateway Timeout")
    )
    result = runner.invoke(app, ["pubchem", "xrefs", "2244"])
    assert result.exit_code == 0
    assert "CHEMBL25" in result.output
    assert "PatentID" in result.output
    assert "(not available)" in result.output
    assert "stops after 30 seconds" in result.output  # the note, on stderr


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


# Hardening audit D2: the output was opened (and emptied) before the input was read,
# so writing results back over the input destroyed it.


@pytest.mark.parametrize(
    ("name", "content"),
    [
        ("compounds.csv", "id,smiles\naspirin,CC(=O)Oc1ccccc1C(=O)O\n"),
        ("compounds.smi", "CCO ethanol\nc1ccccc1 benzene\n"),
    ],
)
def test_batch_refuses_to_write_over_its_input(tmp_path, name, content):
    src = tmp_path / name
    src.write_text(content)
    same = tmp_path / "sub" / ".." / name  # the same file, spelled differently
    (tmp_path / "sub").mkdir()
    for output in (src, same):
        result = runner.invoke(app, ["props", "-i", str(src), "-o", str(output)])
        assert result.exit_code == 2
        assert "is the input file" in result.output
    assert src.read_text() == content


def test_batch_checks_the_input_before_emptying_the_output(tmp_path):
    out = tmp_path / "out.csv"
    out.write_text("precious,data\n")
    bad = tmp_path / "bad.csv"
    bad.write_text("id,structure\n1,CCO\n")
    result = runner.invoke(app, ["props", "-i", str(bad), "-o", str(out)])
    assert result.exit_code == 1
    assert "No SMILES column" in result.output
    assert out.read_text() == "precious,data\n"


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
