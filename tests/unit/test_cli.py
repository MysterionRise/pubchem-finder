"""Tests for the typer CLI."""

import contextlib
import csv
import json
import os
from pathlib import Path

import pytest
from typer.testing import CliRunner

from pubchem_finder import __version__, hints, service
from pubchem_finder.cli import app
from pubchem_finder.models import RemoteHit, RemoteSearchResults

runner = CliRunner()


def test_pubchem_exact_labels_score_column(monkeypatch):
    hit = RemoteHit(
        cid=76973140,
        smiles="CC(=O)Oc1ccccc1C(=O)O",
        title="Acetylsalicylic Acid-d7",
        ecfp4_similarity=0.405,
    )
    stub = RemoteSearchResults(
        query_smiles="CC(=O)Oc1ccccc1C(=O)O",
        search="identity",
        max_records=25,
        truncated=False,
        hits=[hit],
        note=service.IDENTITY_NOTE.format(identity_type="same_connectivity"),
    )
    monkeypatch.setattr(service, "pubchem_exact_match", lambda *args, **kwargs: stub)
    args = [
        "pubchem",
        "exact",
        "CC(=O)Oc1ccccc1C(=O)O",
        "--identity",
        "same_connectivity",
    ]
    output = runner.invoke(app, args).output
    before_first_hit = output.partition("CID 76973140")[0]
    assert "ECFP4" in before_first_hit, output
    assert "similarity" in before_first_hit, output
    assert "below 1" in output, output
    assert "fingerprint tells isotopes apart" in output, output
    assert "Stereoisomers score 1" in output, output


def test_pubchem_exact_empty_results_omit_score_column_legend(monkeypatch):
    stub = RemoteSearchResults(
        query_smiles="CC(=O)Oc1ccccc1C(=O)O",
        search="identity",
        max_records=25,
        truncated=False,
        hits=[],
        note=service.IDENTITY_NOTE.format(identity_type="same_connectivity"),
    )
    monkeypatch.setattr(service, "pubchem_exact_match", lambda *args, **kwargs: stub)
    args = [
        "pubchem",
        "exact",
        "CC(=O)Oc1ccccc1C(=O)O",
        "--identity",
        "same_connectivity",
    ]

    output = runner.invoke(app, args).output

    assert "0 PubChem identity hit(s)" in output, output
    assert "score column: ECFP4 similarity" not in output, output


def test_pubchem_substructure_labels_score_column(monkeypatch):
    hit = RemoteHit(
        cid=1,
        smiles="c1ccc2[nH]ccc2c1",
        title="Indole",
        ecfp4_similarity=1.0,
    )
    stub = RemoteSearchResults(
        query_smiles="c1ccc2[nH]ccc2c1",
        search="substructure",
        max_records=50,
        truncated=False,
        hits=[hit],
        note="",
    )
    monkeypatch.setattr(
        service, "pubchem_substructure_search", lambda *args, **kwargs: stub
    )
    output = runner.invoke(app, ["pubchem", "substructure", "c1ccc2[nH]ccc2c1"]).output
    before_first_hit = output.partition("CID 1 ")[0]
    assert "ECFP4" in before_first_hit, output
    assert "similarity" in before_first_hit, output


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
    assert result.stdout.strip() == "CC(=O)O"
    assert "removed: [Na+]" in result.stderr


def test_standardize_says_when_the_input_is_returned_as_given():
    result = runner.invoke(app, ["standardize", "N.N.Cl[Pt]Cl"])
    assert result.exit_code == 0
    assert "Pt" in result.stdout
    assert "note: Contains Pt" in result.stderr


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
    # `lookup 2244` and `lookup <SMILES>` end here too: say what --by is for.
    assert "--by cid, --by smiles or --by inchikey" in result.output
    # By CID there is nothing to suggest. PubChem answers a CID it does not have
    # with a row that holds the CID alone; without this route the command died on
    # the unmocked request and the assertion below held for an empty output.
    from tests.conftest import property_table

    pubchem_api.get(path__regex=r"/compound/cid/999999999/property/.*").respond(
        200, json=property_table({"CID": 999999999})
    )
    by_cid = runner.invoke(app, ["lookup", "999999999", "--by", "cid"])
    assert by_cid.exit_code == 1
    assert "error: No PubChem compound found for cid='999999999'" in by_cid.output
    assert "--by cid" not in by_cid.output


def test_index_build_requires_one_source(tmp_path):
    result = runner.invoke(app, ["index", "build", "--home", str(tmp_path)])
    assert result.exit_code == 2
    assert "exactly one" in result.output


def test_index_build_rejects_sdf_without_a_chunk_name(tmp_path):
    path = tmp_path / ".sdf"
    path.touch()
    result = runner.invoke(
        app, ["index", "build", "--sdf", str(path), "--home", str(tmp_path / "home")]
    )
    assert result.exit_code == 1
    assert "no chunk name" in result.output


def test_search_reports_invalid_index_mode(built_home, monkeypatch):
    # The mode is read when an index is first opened, and built_home is shared by
    # the whole session: start without the indexes earlier tests opened, or this
    # one passes or fails with the order the files run in.
    monkeypatch.setattr(service, "_indexes", {})
    monkeypatch.setenv("PUBCHEM_FINDER_INDEX_MODE", "bogus")
    result = runner.invoke(app, ["search", "CCO", "--home", str(built_home)])
    assert result.exit_code == 1, result.output
    assert "PUBCHEM_FINDER_INDEX_MODE" in result.output
    assert "auto, memory, or disk" in result.output


def test_search_header_says_when_the_list_was_cut_at_top(built_home):
    """Issue #191: "20 hit(s) >= 0.6" was printed where 44 compounds reach 0.6."""

    def search(top, *more):
        args = ["search", ASPIRIN, "--threshold", "0.1", "--top", top, *more]
        return runner.invoke(app, [*args, "--home", str(built_home)])

    cut = search("2")
    assert cut.exit_code == 0, cut.output
    assert (
        "top 2 of more than 2 hit(s) >= 0.1 tanimoto among 11 compounds" in cut.stdout
    )
    assert "More compounds reach the threshold than the 2 returned" in cut.stderr
    whole = search("50")
    assert whole.exit_code == 0, whole.output
    n = len([line for line in whole.stdout.splitlines() if " CID " in line])
    assert n > 2
    assert whole.stdout.startswith(f"{n} hit(s) >= 0.1 tanimoto")
    assert "note:" not in whole.stderr
    assert json.loads(search("2", "--json").stdout)["truncated"] is True


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
    assert "None of the 1 chunk(s)" in result.output
    assert "holds a compound" in result.output


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
    import shutil

    home = tmp_path / "home"
    assert _sdf_build(home, same_named[0]).exit_code == 0
    result = _sdf_build(home, same_named[1])
    assert result.exit_code == 1
    assert "already holds" in result.output
    assert 2244 in _cids(home)  # vendor A's compounds are still there

    # The way to replace it, as the message gives it: the advice used to stop at
    # "copy the new file over it", after which the same command was refused again.
    first = same_named[0].resolve()
    assert f"copy the new file over it and build from {first}" in result.output
    shutil.copyfile(same_named[1], same_named[0])
    assert _sdf_build(home, same_named[1]).exit_code == 1
    assert _sdf_build(home, same_named[0]).exit_code == 0
    cids = _cids(home)
    assert 1983 in cids
    assert 2244 not in cids


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


def test_index_build_prefers_legacy_spelling_for_unrecorded_chunks(
    tmp_path, sample_sdfs
):
    home = tmp_path / "home"
    modern = tmp_path / "a" / "LIB.sdf.gz"
    legacy = tmp_path / "b" / "LIB.SDF.GZ"
    for path, sample in zip((modern, legacy), sample_sdfs, strict=False):
        path.parent.mkdir()
        path.write_bytes(sample.read_bytes())
    from pubchem_finder.index.builder import SdfSource, build_chunks
    from pubchem_finder.index.layout import IndexPaths

    build_chunks(
        IndexPaths(home),
        [
            SdfSource(name="LIB", local_path=modern),
            SdfSource(name=legacy.name, local_path=legacy),
        ],
        workers=1,
    )
    modern_stats = _as_2_2_0_wrote_it(home, "LIB")
    _as_2_2_0_wrote_it(home, legacy.name)

    result = _sdf_build(home, legacy)

    assert result.exit_code == 0, result.output
    assert (
        f"rebuilding {legacy.name}: its chunk does not record the file" in result.output
    )
    assert json.loads(modern_stats.read_text()).get("source") is None
    assert {1983, 2244} <= _cids(home)


def test_index_build_rebuilds_a_chunk_that_records_no_file(
    tmp_path, same_named, sample_sdfs
):
    """A chunk from 2.2.0 or earlier does not say which file it holds. A draft of
    2.2.2 took a same-named file for that file whenever it was not newer than the
    chunk, and downloads and copies keep the publisher's time: a replaced library
    was recorded as built and never indexed."""
    import shutil

    home = tmp_path / "home"
    library = same_named[0]
    assert _sdf_build(home, library).exit_code == 0
    stats = _as_2_2_0_wrote_it(home)
    shutil.copyfile(sample_sdfs[1], library)  # another version of the library...
    earlier = stats.stat().st_mtime_ns - 5_000_000_000
    os.utime(library, ns=(earlier, earlier))  # ...that looks older than the chunk
    result = _sdf_build(home, library)
    assert result.exit_code == 0, result.output
    assert "rebuilding library: its chunk does not record the file" in result.output
    cids = _cids(home)
    assert 1983 in cids
    assert 2244 not in cids
    assert json.loads(stats.read_text())["source"]["path"] == str(library.resolve())

    again = _sdf_build(home, library)
    assert "index up to date" in again.output
    assert "merging" not in again.output


def _snapshot(folder):
    return {
        str(p.relative_to(folder)): p.read_bytes()
        for p in sorted(folder.rglob("*"))
        if p.is_file()
    }


@pytest.mark.parametrize("content", ["empty", "error page", "text ids"])
@pytest.mark.parametrize("with_another_chunk", [False, True])
def test_file_without_usable_records_does_not_replace_its_chunk(
    tmp_path, sample_sdfs, content, with_another_chunk
):
    """A failed export (0 bytes, a saved error page) or one whose IDs changed form
    was fingerprinted into an empty chunk that replaced the good one. With one
    library the index then stayed "ready" on fingerprints without records, and
    every later build failed; with two, the library left the index silently."""
    import shutil

    from rdkit import Chem

    from tests.conftest import CHUNK_A, write_sdf

    home = tmp_path / "home"
    library = write_sdf(tmp_path / "library.sdf", CHUNK_A, gz=False)
    files = [library]
    if with_another_chunk:
        files.append(shutil.copyfile(sample_sdfs[1], tmp_path / "other.sdf.gz"))
    assert _sdf_build(home, *files).exit_code == 0
    before = _snapshot(home / "chunks")

    if content == "empty":
        library.write_bytes(b"")
    elif content == "error page":
        library.write_text("<html><body>404 Not Found</body></html>\n")
    else:
        with Chem.SDWriter(str(library)) as writer:
            for i, smiles in enumerate(["CCO", "CCN", "CCC"], 1):
                mol = Chem.MolFromSmiles(smiles)
                mol.SetProp("_Name", f"ABC-{i:04d}")
                writer.write(mol)
    result = _sdf_build(home, *files)
    assert result.exit_code == 1, result.output
    assert "has no record the index can hold" in result.output
    assert f"Chunk 'library' keeps the {len(CHUNK_A)} compounds it has" in result.output
    if content == "text ids":
        assert "skipped 'ABC-0001' and any others like it" in result.output
    assert _snapshot(home / "chunks") == before  # the chunk, and nothing left over
    found = runner.invoke(
        app, ["search", "CC(=O)Oc1ccccc1C(=O)O", "--home", str(home), "--json"]
    )
    hit = json.loads(found.output)["hits"][0]
    assert (hit["cid"], hit["smiles"]) == (2244, "CC(=O)Oc1ccccc1C(=O)O")


@pytest.mark.parametrize("with_another_chunk", [False, True])
def test_the_way_out_of_an_empty_file_is_one_that_works(
    tmp_path, sample_sdfs, with_another_chunk
):
    """The refusal said "delete library.json and build again" whatever the index
    held. With one library that build fails too ("nothing to finalize"), and the
    index went on answering from the old fingerprints, without a way forward."""
    import shutil

    from tests.conftest import CHUNK_A, write_sdf

    home = tmp_path / "home"
    library = write_sdf(tmp_path / "library.sdf", CHUNK_A, gz=False)
    files = [library]
    if with_another_chunk:
        files.append(shutil.copyfile(sample_sdfs[1], tmp_path / "other.sdf.gz"))
    assert _sdf_build(home, *files).exit_code == 0
    library.write_bytes(b"")
    refused = _sdf_build(home, *files)
    assert refused.exit_code == 1
    stats = home / "chunks" / "library.json"

    if with_another_chunk:
        assert f"delete {stats} and build again" in refused.output
        stats.unlink()
        result = _sdf_build(home, *files)
        assert result.exit_code == 0, result.output
        # The chunk is complete, and empty: "no longer complete ... Build from
        # its file again to put it back" was the wrong thing to say about it.
        assert "warning: library: in the index, but its file has no record" in (
            result.output
        )
        assert "no longer complete" not in result.output
        cids = _cids(home)
        assert 1983 in cids
        assert 2244 not in cids
        return
    advice = f"delete {home / 'index'} and {home / 'chunks'}"
    assert f"It is all this index holds: to remove the index instead, {advice}" in (
        refused.output
    )
    assert "build again" not in refused.output
    # Deleting the stats file all the same: the build says what state that leaves.
    stats.unlink()
    result = _sdf_build(home, *files)
    assert result.exit_code == 1
    assert "nothing to finalize" in result.output
    assert "The index still answers from the earlier build of library" in result.output
    assert f"or {advice} to remove the index" in result.output
    shutil.rmtree(home / "index")
    shutil.rmtree(home / "chunks")
    status = runner.invoke(app, ["index", "status", "--home", str(home), "--json"])
    assert json.loads(status.output)["state"] == "not_built"


def test_an_empty_chunk_may_be_replaced_by_another_empty_one(tmp_path):
    """The refusal is for a chunk that holds compounds; a file that never gave any
    is simply read again."""
    from rdkit import Chem

    home = tmp_path / "home"
    lib = tmp_path / "chembl.sdf"
    for name in ("CHEMBL25", "CHEMBL26"):
        mol = Chem.MolFromSmiles("CCO")
        mol.SetProp("_Name", name)
        with Chem.SDWriter(str(lib)) as writer:
            writer.write(mol)
        result = _sdf_build(home, lib)
        assert result.exit_code == 1
        assert "None of the 1 chunk(s)" in result.output
        assert f"skipped {name!r}" in result.output
        assert "keeps the" not in result.output


def test_rebuild_from_a_damaged_file_leaves_the_chunk_as_it_was(
    tmp_path, same_named, sample_sdfs
):
    """The real write_chunk failing part-way through a file: 2.2.1 had deleted the
    old chunk by then, and its tests only ever failed before the first record."""
    home = tmp_path / "home"
    library = same_named[0]
    assert _sdf_build(home, library).exit_code == 0
    before = _snapshot(home / "chunks")

    data = sample_sdfs[1].read_bytes()
    library.write_bytes(data[: len(data) // 2])  # a download that stopped half-way
    result = _sdf_build(home, library)
    assert result.exit_code == 1
    assert f"{library.resolve()} is damaged or incomplete" in result.output
    assert _snapshot(home / "chunks") == before  # no .tmp files either
    assert 2244 in _cids(home)


def test_moved_and_changed_file_says_how_to_replace_its_chunk(
    tmp_path, same_named, sample_sdfs
):
    """The advice was "Rename one of them" although only one file was left, and
    renaming it indexes the library a second time next to its old version."""
    import shutil

    home = tmp_path / "home"
    library = same_named[0]
    assert _sdf_build(home, library).exit_code == 0
    recorded = library.resolve()
    moved = tmp_path / "elsewhere" / library.name
    moved.parent.mkdir()
    shutil.copyfile(sample_sdfs[1], moved)
    library.unlink()
    result = _sdf_build(home, moved)
    assert result.exit_code == 1
    assert f"was built from {recorded}, which is no longer there" in result.output
    assert f"put the file at {recorded} and build again" in result.output
    assert "Rename one of them" not in result.output
    assert 2244 in _cids(home)


def test_merge_says_when_it_leaves_out_a_chunk_that_is_no_longer_complete(
    tmp_path, same_named, sample_sdfs
):
    """A rebuild killed inside its swap leaves the chunk without its stats file; a
    later build of something else merged the library out without a word."""
    import shutil

    home = tmp_path / "home"
    library = same_named[0]
    other = shutil.copyfile(sample_sdfs[1], tmp_path / "other.sdf.gz")
    assert _sdf_build(home, library, other).exit_code == 0
    (home / "chunks" / "library.json").unlink()
    result = _sdf_build(home, other)
    assert result.exit_code == 0, result.output
    assert "warning: library: in the index, but no longer complete" in result.output
    assert 2244 not in _cids(home)


def _built_chunks(home):
    from pubchem_finder.index.layout import IndexPaths

    return IndexPaths(home).completed_chunks()


def test_chunks_that_differ_only_by_case_stay_buildable(tmp_path, sample_sdfs):
    """An index made on a case-sensitive filesystem before 2.2.2 may hold chunks
    Library and library. A draft of 2.2.2 kept one of them per folded name and
    refused the other's own file with "Rename the file"."""
    import shutil

    from pubchem_finder.index.builder import SdfSource, _check_unique_names

    upper = tmp_path / "a" / "Library.sdf.gz"
    lower = tmp_path / "b" / "library.sdf.gz"
    stats = []
    for path, src in ((upper, sample_sdfs[0]), (lower, sample_sdfs[1])):
        path.parent.mkdir()
        shutil.copyfile(src, path)
        home = tmp_path / f"home-{path.parent.name}"
        assert _sdf_build(home, path).exit_code == 0
        stats += _built_chunks(home)
    assert sorted(c.name for c in stats) == ["Library", "library"]

    sources = [
        SdfSource(name="Library", local_path=upper),
        SdfSource(name="library", local_path=lower),
    ]
    for given in ([sources[0]], [sources[1]], sources):
        assert _check_unique_names(given, stats) == given


def test_a_chunk_without_a_recorded_file_knows_it_under_another_case(
    tmp_path, sample_sdfs
):
    """On macOS and Windows `library.sdf` and `Library.sdf` are one file. A chunk
    from 2.2.0 or earlier records no file to compare with, so the file typed in
    another case was refused with "Rename the file"; 2.2.1 rebuilt it."""
    import shutil
    from dataclasses import replace

    from pubchem_finder.index.builder import SdfSource, _check_unique_names

    typed = shutil.copyfile(sample_sdfs[0], tmp_path / "LIBRARY.sdf.gz")
    # Where names keep their case, a hard link under the chunk's spelling; where
    # they don't, the file answers to that spelling already.
    with contextlib.suppress(FileExistsError):
        os.link(typed, tmp_path / "library.sdf.gz")
    home = tmp_path / "home"
    assert _sdf_build(home, tmp_path / "library.sdf.gz").exit_code == 0
    _as_2_2_0_wrote_it(home)
    src = SdfSource(name="LIBRARY", local_path=typed)
    assert _check_unique_names([src], _built_chunks(home)) == [
        replace(src, name="library")
    ]


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


def test_a_rebuild_on_an_index_merged_by_2_2_1_gives_its_own_reason(
    tmp_path, same_named, sample_sdfs
):
    """ "merging once more: ... 2.2.1 or earlier" was printed under "rebuilding
    library: ... changed since it was indexed", as a second reason for one merge."""
    import shutil

    home = tmp_path / "home"
    library = same_named[0]
    assert _sdf_build(home, library).exit_code == 0
    for path, field in (
        (home / "index" / "manifest.json", "chunk_builds"),
        (home / "chunks" / "library.json", "build_id"),
    ):  # as 2.2.1 wrote them
        data = json.loads(path.read_text())
        del data[field]
        path.write_text(json.dumps(data))
    shutil.copyfile(sample_sdfs[1], library)
    result = _sdf_build(home, library)
    assert result.exit_code == 0, result.output
    assert "rebuilding library: " in result.output
    assert "merging once more" not in result.output
    assert 1983 in _cids(home)


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
    assert "merging once more: this index was merged by pubchem-finder 2.2.1" in (
        result.output
    )
    assert "rebuilding" not in result.output
    cids = _cids(home)
    assert 1983 in cids
    assert 2244 not in cids
    # The chunk gets a build id as it is merged, so that a rebuild by 2.2.1 (which
    # writes none) is not taken for the build the index holds.
    stats_path = home / "chunks" / "library.json"
    stats = json.loads(stats_path.read_text())
    manifest = json.loads((home / "index" / "manifest.json").read_text())
    assert stats["build_id"]
    assert manifest["chunk_builds"] == {"library": stats["build_id"]}
    again = _sdf_build(home, library)
    assert "index up to date" in again.output
    assert "merging" not in again.output

    del stats["build_id"]  # 2.2.1 rebuilds the chunk and fails before its merge
    stats_path.write_text(json.dumps(stats))
    assert "merging 1 chunk(s)" in _sdf_build(home, library).output


# Found reviewing the repair above, on real files: a rebuilt chunk's records go live
# as soon as its files are swapped in, while the merged index holds the fingerprints
# of the file it replaced until the merge is done. When the build failed or was
# stopped in between, a search paired the two: with IDs that mean another compound in
# the new export, aspirin's fingerprint came back with caffeine's structure at 1.000,
# and `index status` said ready.

ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"
CAFFEINE = "Cn1cnc2c1c(=O)n(C)c(=O)n2C"
PARACETAMOL = "CC(=O)Nc1ccc(O)cc1"
EXPORT_1 = [(1, ASPIRIN, "aspirin", "new"), (2, "c1ccccc1", "benzene", "new")]
EXPORT_2 = [(1, CAFFEINE, "caffeine", "new"), (5, PARACETAMOL, "paracetamol", "new")]


def _search(home, smiles, *more):
    found = runner.invoke(
        app, ["search", smiles, "--home", str(home), "--threshold", "0.9", *more]
    )
    assert found.exit_code == 0, found.output
    return found


def _status(home):
    result = runner.invoke(app, ["index", "status", "--home", str(home), "--json"])
    return json.loads(result.output)


@pytest.fixture
def renumbered(tmp_path):
    """An index of export 1, with the file then replaced by export 2."""
    from tests.conftest import write_sdf

    home = tmp_path / "home"
    library = write_sdf(tmp_path / "library.sdf.gz", EXPORT_1)
    args = ["index", "build", "--home", str(home), "--workers", "1"]
    result = runner.invoke(app, [*args, "--sdf", str(library), "--with-substructure"])
    assert result.exit_code == 0, result.output
    write_sdf(library, EXPORT_2)
    return home, library


def test_rebuilt_chunk_is_not_read_with_the_fingerprints_it_replaced(
    renumbered, monkeypatch
):
    from pubchem_finder.index import builder

    home, library = renumbered
    real = builder.sort_fpsim2_file
    monkeypatch.setattr(
        builder,
        "sort_fpsim2_file",
        lambda *a, **k: (_ for _ in ()).throw(OSError(28, "No space left on device")),
    )
    failed = _sdf_build(home, library)
    assert failed.exit_code == 1
    assert "built library: 2 compounds" in failed.output  # the chunk was replaced
    assert "No space left on device" in failed.output
    assert "The current index is unchanged" not in failed.output
    assert "still holds the earlier build of library" in failed.output

    # The index answers from export 1's fingerprints. ID 1 is caffeine now, so the
    # hit must not carry a structure, and the result has to say why.
    hits = json.loads(_search(home, ASPIRIN, "--json").output)
    assert [(h["cid"], h["smiles"]) for h in hits["hits"]] == [(1, None)]
    assert "library was indexed again after the last merge" in hits["warning"]
    # Quoted as the message quotes it: a temporary directory may have a space.
    quoted = hints._quote(str(library.resolve()))
    assert f"index build --sdf {quoted}" in hits["warning"]
    text = _search(home, ASPIRIN)
    assert "warning: The index is behind its files" in text.output
    assert CAFFEINE not in text.output

    sub = runner.invoke(
        app, ["substructure", "c1ccccc1", "--home", str(home), "--json"]
    )
    found = json.loads(sub.output)
    assert found["hits"] == []  # not benzene's ID with another structure
    assert "library was indexed again" in found["warning"]
    assert "its compounds cannot be checked and are left out" in found["warning"]
    # In words: raising --max would not help, the warning says what will.
    text = runner.invoke(app, ["substructure", "c1ccccc1", "--home", str(home)])
    assert "0 hit(s) (more may exist; see the warning); 2 of 2" in text.output
    assert "raise --max" not in text.output
    # Two candidates passed the screen and neither could be checked: the answer
    # used to call itself exhaustive.
    assert (found["n_candidates"], found["exhaustive"]) == (2, False)

    status = _status(home)
    assert status["state"] == "ready"
    assert status["unmerged_chunks"] == ["library"]
    assert "The index is behind its files" in status["hint"]
    # The text form skipped every list, and so the one line that says which.
    text = runner.invoke(app, ["index", "status", "--home", str(home)])
    assert "unmerged_chunks      library" in text.output

    # Batch mode has no column for the warning: it wrote the hit without its
    # structure, "0 error(s)", and nothing else.
    queries = home.parent / "queries.csv"
    queries.write_text(f"smiles,name\n{ASPIRIN},a\n{ASPIRIN},b\n")
    for command, query in (("search", ASPIRIN), ("substructure", "c1ccccc1")):
        queries.write_text(f"smiles,name\n{query},a\n{query},b\n")
        args = [command, "--input", str(queries), "--output", str(home.parent / "out")]
        batch = runner.invoke(app, [*args, "--home", str(home)])
        assert batch.exit_code == 0, batch.output
        assert batch.output.count("warning: The index is behind its files") == 1

    monkeypatch.setattr(builder, "sort_fpsim2_file", real)
    again = _sdf_build(home, library)
    assert again.exit_code == 0, again.output
    assert _status(home)["unmerged_chunks"] == []
    hits = json.loads(_search(home, CAFFEINE, "--json").output)
    assert [(h["cid"], h["smiles"]) for h in hits["hits"]] == [(1, CAFFEINE)]
    assert hits["warning"] is None
    assert json.loads(_search(home, ASPIRIN, "--json").output)["hits"] == []


def test_a_running_server_sees_a_chunk_change_under_it(renumbered):
    """One LocalIndex for the life of the process, as the MCP server keeps it."""
    from pubchem_finder import service
    from pubchem_finder.index.builder import SdfSource, build_chunks, finalize
    from pubchem_finder.index.layout import IndexPaths

    home, library = renumbered
    paths = IndexPaths(home)

    def top(smiles):
        result = service.similarity_search(smiles, threshold=0.9, home=home)
        return [(h.cid, h.smiles) for h in result.hits], result.warning

    assert top(ASPIRIN) == ([(1, ASPIRIN)], None)
    build_chunks(paths, [SdfSource(name="library", local_path=library)])  # no merge
    hits, warning = top(ASPIRIN)
    assert hits == [(1, None)]
    assert warning is not None
    assert service.index_status(home).unmerged_chunks == ["library"]

    finalize(paths)
    assert top(CAFFEINE) == ([(1, CAFFEINE)], None)
    assert top(ASPIRIN) == ([], None)
    assert service.index_status(home).unmerged_chunks == []


def test_a_merge_that_completes_under_a_search_does_not_lend_it_structures(
    renumbered, monkeypatch
):
    """The hits were looked up against the manifest on disk at that moment. When a
    build finished between the search and the lookup, aspirin's hit from the old
    fingerprints came back as caffeine at 1.000, with no warning: once in 763
    searches run next to a rebuild loop."""
    from rdkit import Chem

    from pubchem_finder import service
    from pubchem_finder.index.builder import SdfSource, build_chunks, finalize
    from pubchem_finder.index.engine import LocalIndex
    from pubchem_finder.index.layout import IndexPaths

    home, library = renumbered
    paths = IndexPaths(home)
    index = service.get_index(home)

    def rebuild():
        build_chunks(paths, [SdfSource(name="library", local_path=library)])
        finalize(paths)

    found = index.search_merged(Chem.MolFromSmiles(ASPIRIN), 0.9, 5)
    pairs, merged = found.pairs, found.merged
    assert pairs == [(1, 1.0)]
    rebuild()
    assert index.metadata([1])[1]["smiles"] == CAFFEINE  # what ID 1 is now
    assert index.metadata([1], merged) == {}  # and not what the hit was
    assert list(index.unmerged(merged)) == ["library"]
    assert index.unmerged() == {}

    # The same through the service, for both searches. The index holds export 2
    # now (ID 1 is caffeine); each search is overtaken by a build of the other
    # export, in which ID 1 has a carbonyl too.
    from tests.conftest import write_sdf

    def overtaken(method, run):
        real = getattr(LocalIndex, method)

        def then_a_build_finishes(self, *args, **kwargs):
            found = real(self, *args, **kwargs)
            rebuild()
            return found

        with monkeypatch.context() as patched:
            patched.setattr(LocalIndex, method, then_a_build_finishes)
            return run()

    write_sdf(library, EXPORT_1)
    similar = overtaken(
        "search_merged", lambda: service.similarity_search(CAFFEINE, 0.9, home=home)
    )
    assert [(h.cid, h.smiles) for h in similar.hits] == [(1, None)]
    assert "library was indexed again" in similar.warning

    write_sdf(library, EXPORT_2)
    containing = overtaken(
        "_screen", lambda: service.substructure_search("C=O", home=home)
    )
    assert containing.n_candidates == 1  # aspirin, in the fingerprints searched
    assert containing.hits == []  # not caffeine, which has its ID now
    assert containing.exhaustive is False
    assert "library was indexed again" in containing.warning

    after = service.similarity_search(CAFFEINE, 0.9, home=home)
    assert ([(h.cid, h.smiles) for h in after.hits], after.warning) == (
        [(1, CAFFEINE)],
        None,
    )


def test_a_lookup_never_returns_rows_it_could_not_verify(renumbered, monkeypatch):
    """The lookup read again when a chunk changed under it, three times at most,
    and then returned the last read as it was. With three or more files being
    rebuilt (or one rebuilt three times) during one lookup, a 1.000 hit carried
    another build's structure: 12 times in 43,349 searches next to a real build
    of 24 files. A row now counts only if its own chunk was the merged build both
    before and after it was read."""
    from pubchem_finder import service
    from pubchem_finder.index.builder import SdfSource, build_chunks
    from pubchem_finder.index.engine import LocalIndex
    from pubchem_finder.index.layout import IndexPaths
    from tests.conftest import CHUNK_B, write_sdf

    home, library = renumbered
    paths = IndexPaths(home)
    other = write_sdf(home.parent / "other.sdf.gz", CHUNK_B)
    write_sdf(library, EXPORT_1)
    assert _sdf_build(home, library, other).exit_code == 0
    exports = [EXPORT_2, EXPORT_1]
    real = LocalIndex._metadata
    reads = []

    def a_build_replaces_the_library_during_every_read(self, cids, chunks):
        found = real(self, cids, chunks)
        write_sdf(library, exports[len(reads) % 2])
        build_chunks(paths, [SdfSource(name="library", local_path=library)])
        reads.append(chunks)
        return found

    monkeypatch.setattr(
        LocalIndex, "_metadata", a_build_replaces_the_library_during_every_read
    )
    result = service.similarity_search(ASPIRIN, threshold=0.9, home=home)
    assert [(h.cid, h.smiles) for h in result.hits] == [(1, None)]
    assert "library was indexed again" in result.warning
    assert len(reads) == 1  # nothing to gain from reading again
    # The chunk that did not change is read as before, in the same lookup.
    index = service.get_index(home)
    rows = index.metadata([1, 2519], index.merged)
    assert list(rows) == [2519]
    assert rows[2519]["smiles"] == CAFFEINE


def test_a_file_name_duckdb_reads_as_a_pattern_is_read_as_a_name(tmp_path):
    """DuckDB takes [ ] * and ? in a path for a pattern: asked for
    "lib[1].parquet" it opened "lib1.parquet". A hit from one library came back
    without its structure, or with the one the other library has under that ID."""
    from pubchem_finder.index.layout import chunk_of, duckdb_path
    from tests.conftest import write_sdf

    assert duckdb_path("/data/a*b?c[1].parquet") == "/data/a[*]b[?]c[[]1].parquet"
    assert duckdb_path("/data/plain-name_1.parquet") == "/data/plain-name_1.parquet"
    if os.sep == "/":  # no pattern stands for a backslash: such a path is left alone
        assert duckdb_path("/data/a\\b[1].parquet") == "/data/a\\b[1].parquet"
    # Back from the path of a row to its chunk, also for a file called ".sdf".
    assert chunk_of("/data/chunks/lib.v2.parquet") == "lib.v2"
    assert chunk_of("/data/chunks/.parquet") == ""
    home = tmp_path / "home [x]"
    first = write_sdf(tmp_path / "lib[1].sdf", [(11, ASPIRIN, "a", "new")], gz=False)
    second = write_sdf(tmp_path / "lib1.sdf", [(12, CAFFEINE, "c", "new")], gz=False)
    assert _sdf_build(home, first, second).exit_code == 0
    hits = json.loads(_search(home, ASPIRIN, "--json").output)["hits"]
    assert [(h["cid"], h["smiles"]) for h in hits] == [(11, ASPIRIN)]
    hits = json.loads(_search(home, CAFFEINE, "--json").output)["hits"]
    assert [(h["cid"], h["smiles"]) for h in hits] == [(12, CAFFEINE)]
    # The check of IDs reads the right files too: with ID 11 in both, the build
    # compared lib1 with itself and merged aspirin and caffeine under one ID.
    write_sdf(second, [(11, CAFFEINE, "c", "new")], gz=False)
    result = _sdf_build(home, first, second)
    assert result.exit_code == 1
    assert "ID 11 stands for different compounds" in result.output
    # Under the names the files have, not as they are spelled for DuckDB.
    assert f"{CAFFEINE} in lib1 and {ASPIRIN} in lib[1]" in result.output
    assert "Give the records of lib1 other IDs" in result.output


@pytest.mark.skipif(os.sep != "/", reason="a backslash separates folders here")
def test_a_home_with_a_backslash_and_a_bracket_in_its_name_is_searched(tmp_path):
    """Escaping the bracket for DuckDB made the backslash an escape too, and a home
    that had worked ended in a traceback."""
    from tests.conftest import write_sdf

    home = tmp_path / "back\\slash [1]" / "home"
    library = write_sdf(tmp_path / "lib.sdf", [(11, ASPIRIN, "a", "new")], gz=False)
    assert _sdf_build(home, library).exit_code == 0
    hits = json.loads(_search(home, ASPIRIN, "--json").output)["hits"]
    assert [(h["cid"], h["smiles"]) for h in hits] == [(11, ASPIRIN)]


def test_a_record_with_a_structure_is_shown_before_one_without(tmp_path):
    """One ID in two files, in one of them as a record without atoms: the hit
    showed whichever row was read first, and so sometimes no structure at all."""
    from rdkit import Chem

    from tests.conftest import write_sdf

    home = tmp_path / "home"
    blank = tmp_path / "a_blank.sdf"
    with Chem.SDWriter(str(blank)) as writer:
        mol = Chem.Mol()
        mol.SetProp("_Name", "11")
        writer.write(mol)
        other = Chem.MolFromSmiles("CCO")
        other.SetProp("_Name", "12")
        writer.write(other)
    library = write_sdf(tmp_path / "z_lib.sdf", [(11, ASPIRIN, "a", "new")], gz=False)
    assert _sdf_build(home, blank, library).exit_code == 0
    hits = json.loads(_search(home, ASPIRIN, "--json").output)["hits"]
    assert [(h["cid"], h["smiles"]) for h in hits] == [(11, ASPIRIN)]


def test_a_chunk_whose_stats_are_gone_is_not_called_indexed_again(renumbered):
    from pubchem_finder import service

    home, _ = renumbered
    # A process that has read the stats before must notice too (the MCP server
    # keeps them between searches).
    assert service.index_status(home).unmerged_chunks == []
    (home / "chunks" / "library.json").unlink()
    found = service.similarity_search(ASPIRIN, threshold=0.9, home=home)
    assert [(h.cid, h.smiles) for h in found.hits] == [(1, None)]
    assert "library is no longer complete in the chunks directory" in found.warning
    assert "indexed again" not in found.warning
    # No file is on record for it, so the hint cannot name one.
    assert "the `index build` command they were built with" in found.warning
    assert service.index_status(home).unmerged_chunks == ["library"]


def test_failed_build_merges_what_it_had_already_replaced(renumbered, tmp_path):
    """One damaged file among several: the others were rebuilt, nothing was merged,
    and the only message was about the damaged file."""
    from tests.conftest import CHUNK_B, write_sdf

    home, library = renumbered
    other = write_sdf(tmp_path / "other.sdf.gz", CHUNK_B)
    data = other.read_bytes()
    other.write_bytes(data[: len(data) // 2])
    failed = _sdf_build(home, library, other)
    assert failed.exit_code == 1
    assert "library: indexed again before the build failed; merging" in failed.output
    assert failed.output.rstrip().splitlines()[-1].startswith("error: ")
    assert "is damaged or incomplete" in failed.output.rstrip().splitlines()[-1]

    assert _status(home)["unmerged_chunks"] == []
    hits = json.loads(_search(home, CAFFEINE, "--json").output)
    assert [(h["cid"], h["smiles"]) for h in hits["hits"]] == [(1, CAFFEINE)]
    assert hits["warning"] is None
    found = json.loads(_search(home, PARACETAMOL, "--json").output)
    assert [h["cid"] for h in found["hits"]] == [5]


def test_failed_build_does_not_merge_out_a_chunk_that_lost_its_stats(
    tmp_path, sample_sdfs
):
    """Only chunks the failed build replaced are merged; one whose stats file is
    gone would be merged out of the index."""
    import shutil

    home = tmp_path / "home"
    library = shutil.copyfile(sample_sdfs[0], tmp_path / "library.sdf.gz")
    assert _sdf_build(home, library).exit_code == 0
    (home / "chunks" / "library.json").unlink()
    broken = tmp_path / "broken.sdf.gz"
    broken.write_bytes(sample_sdfs[1].read_bytes()[:200])
    failed = _sdf_build(home, broken)
    assert failed.exit_code == 1
    assert "merging" not in failed.output
    assert _status(home)["chunks_in_index"] == 1


def test_failed_build_that_replaced_nothing_merges_nothing(
    tmp_path, same_named, sample_sdfs
):
    """New chunks wait for a build that succeeds, as they always did: a download
    that fails at chunk 37 of 359 must not start a merge of the first 36."""
    home = tmp_path / "home"
    assert _sdf_build(home, same_named[0]).exit_code == 0
    new = tmp_path / "new.sdf.gz"
    new.write_bytes(sample_sdfs[1].read_bytes())
    broken = tmp_path / "broken.sdf.gz"
    broken.write_bytes(sample_sdfs[1].read_bytes()[:200])
    failed = _sdf_build(home, new, broken)
    assert failed.exit_code == 1
    assert "merging" not in failed.output
    assert _status(home)["chunks_in_index"] == 1
    assert _status(home)["unmerged_chunks"] == []


# A hit's structure is looked up by its ID alone. Two vendor files that each number
# their records 1..N were merged without a word, and a search for aspirin then
# returned paracetamol's structure at similarity 1.000; a benzene substructure
# search found 2 of the 4 compounds that contain it.

VENDOR_A = [(1, ASPIRIN, "aspirin", "title"), (2, CAFFEINE, "caffeine", "title")]
VENDOR_B = [(1, PARACETAMOL, "paracetamol", "title"), (2, "c1ccccc1", "x", "title")]


def test_files_that_reuse_ids_for_other_compounds_are_not_merged(tmp_path):
    from tests.conftest import write_sdf

    home = tmp_path / "home"
    a = write_sdf(tmp_path / "vendorA.sdf", VENDOR_A, gz=False)
    b = write_sdf(tmp_path / "vendorB.sdf", VENDOR_B, gz=False)
    assert _sdf_build(home, a).exit_code == 0
    result = _sdf_build(home, a, b)
    assert result.exit_code == 1
    message = result.output.rstrip().splitlines()[-1]
    assert "error: ID 1 stands for different compounds" in message
    assert f"{ASPIRIN} in vendorA and {PARACETAMOL} in vendorB" in message
    assert "Give the records of vendorB other IDs" in message  # the new file
    assert str(home / "chunks" / "vendorB.json") in message
    assert "Nothing was merged" in message
    # The index built before is as it was, and right.
    hits = json.loads(_search(home, ASPIRIN, "--json").output)["hits"]
    assert [(h["cid"], h["smiles"]) for h in hits] == [(1, ASPIRIN)]
    assert _status(home)["chunks_in_index"] == 1

    # Renumbered, the second file goes in.
    write_sdf(b, [(101, PARACETAMOL, "paracetamol", "title")], gz=False)
    assert _sdf_build(home, a, b).exit_code == 0
    hits = json.loads(_search(home, PARACETAMOL, "--json").output)["hits"]
    assert [(h["cid"], h["smiles"]) for h in hits] == [(101, PARACETAMOL)]


def test_an_id_used_twice_in_one_file_is_not_merged(tmp_path):
    from tests.conftest import write_sdf

    library = write_sdf(
        tmp_path / "library.sdf",
        [(700001, "CCO", "ethanol", "title"), (700001, "CCN", "ethylamine", "title")],
        gz=False,
    )
    result = _sdf_build(tmp_path / "home", library)
    assert result.exit_code == 1
    assert "ID 700001 stands for different compounds" in result.output
    assert "CCO and CCN, both in library" in result.output


def test_a_library_numbered_from_one_does_not_go_next_to_pubchems_chunks(
    tmp_path, sample_sdfs
):
    """PubChem's CID 702 is ethanol; a library whose record 702 is something else
    would be shown with ethanol's structure, or ethanol with the library's."""
    from pubchem_finder.index.builder import SdfSource, build_chunks, finalize
    from pubchem_finder.index.layout import IndexPaths, chunk_name
    from tests.conftest import write_sdf

    home = tmp_path / "home"
    paths = IndexPaths(home)
    # A chunk as a download leaves it: no recorded file.
    name = chunk_name(sample_sdfs[0])
    build_chunks(paths, [SdfSource(name=name, fetch=lambda: sample_sdfs[0])])
    finalize(paths)
    library = write_sdf(
        tmp_path / "library.sdf", [(702, CAFFEINE, "caffeine", "title")], gz=False
    )
    result = _sdf_build(home, library)
    assert result.exit_code == 1
    assert "ID 702 stands for different compounds" in result.output
    assert "an ID is read as a PubChem CID" in result.output


def test_every_file_is_checked_not_only_the_first(tmp_path):
    from tests.conftest import write_sdf

    a = write_sdf(tmp_path / "a.sdf", [(1, ASPIRIN, "aspirin", "new")], gz=False)
    b = write_sdf(tmp_path / "b.sdf", [(2, CAFFEINE, "caffeine", "new")], gz=False)
    c = write_sdf(tmp_path / "c.sdf", [(2, PARACETAMOL, "x", "new")], gz=False)
    result = _sdf_build(tmp_path / "home", a, b, c)
    assert result.exit_code == 1
    assert "ID 2 stands for different compounds" in result.output
    assert f"{CAFFEINE} in b and {PARACETAMOL} in c" in result.output


def test_an_id_used_twice_in_a_later_file_is_not_merged(tmp_path):
    from tests.conftest import write_sdf

    a = write_sdf(tmp_path / "a.sdf", [(1, ASPIRIN, "aspirin", "title")], gz=False)
    b = write_sdf(
        tmp_path / "b.sdf",
        [(7, "CCO", "ethanol", "title"), (7, "CCN", "ethylamine", "title")],
        gz=False,
    )
    result = _sdf_build(tmp_path / "home", a, b)
    assert result.exit_code == 1
    assert "CCO and CCN, both in b" in result.output


def test_the_advice_names_the_file_that_is_new_or_changed(tmp_path):
    """Whatever the files are called: the advice was to renumber, or delete the
    chunk of, the file that had been in the index all along."""
    from tests.conftest import write_sdf

    home = tmp_path / "home"
    merged = write_sdf(tmp_path / "zeta.sdf", VENDOR_A, gz=False)
    new = write_sdf(tmp_path / "alpha.sdf", VENDOR_B, gz=False)
    assert _sdf_build(home, merged).exit_code == 0
    result = _sdf_build(home, merged, new)
    assert result.exit_code == 1
    assert "Give the records of alpha other IDs" in result.output
    assert str(home / "chunks" / "alpha.json") in result.output

    # Both in the index, and then one of them changes into a clash.
    write_sdf(new, [(101, PARACETAMOL, "paracetamol", "title")], gz=False)
    assert _sdf_build(home, merged, new).exit_code == 0
    write_sdf(merged, [(101, CAFFEINE, "caffeine", "title")], gz=False)
    result = _sdf_build(home, merged, new)
    assert result.exit_code == 1
    assert "ID 101 stands for different compounds" in result.output
    assert "Give the records of zeta other IDs" in result.output
    assert str(home / "chunks" / "zeta.json") in result.output


def test_every_shared_id_is_checked_however_many(tmp_path, monkeypatch):
    """Only the first 20,000 shared IDs were compared: two exports that agree up
    to there and differ after it were merged, with nothing said about a limit."""
    from rdkit import Chem

    from pubchem_finder.index import builder

    monkeypatch.setattr(builder, "ID_CHECK_BATCH", 2)
    acids = ["CC(=O)O", "CCC(=O)O", "CCCC(=O)O", "CCCCC(=O)O", "OC(=O)c1ccccc1"]

    def export(path, smiles):  # no InChIKey tag: only RDKit can compare these
        with Chem.SDWriter(str(path)) as writer:
            for number, text in enumerate(smiles, start=1):
                mol = Chem.MolFromSmiles(text)
                mol.SetProp("_Name", str(number))
                writer.write(mol)
        return path

    first = export(tmp_path / "first.sdf", acids)
    # The same acids as their anions, except the last record.
    anions = [text.replace("(=O)O", "(=O)[O-]") for text in acids[:4]]
    second = export(tmp_path / "second.sdf", [*anions, "CCN"])
    result = _sdf_build(tmp_path / "home", first, second)
    assert result.exit_code == 1
    assert "ID 5 stands for different compounds" in result.output
    assert "O=C(O)c1ccccc1 in first and CCN in second" in result.output
    export(second, [*anions, "[O-]C(=O)c1ccccc1"])
    monkeypatch.setattr(builder, "ID_CHECK_ANNOUNCE", 2)
    merged = _sdf_build(tmp_path / "home", first, second)
    assert merged.exit_code == 0
    # A long comparison says that it is running, and how far it is.
    assert "comparing 5 records of first that share an ID" in merged.output
    assert "compared 4 of 5" in merged.output

    # What one merge has compared, the next does not compare again: adding a
    # file of two compounds to an index of 300,000 shared IDs took 88 s each time.
    asked = []
    real = builder._same_compound
    monkeypatch.setattr(
        builder, "_same_compound", lambda a, b: asked.append((a, b)) or real(a, b)
    )
    third = export(tmp_path / "third.sdf", ["CCCCCCCC"] * 5)
    result = _sdf_build(tmp_path / "home", first, second, third)
    assert result.exit_code == 1  # ID 1 is octane there: compared, and refused
    assert asked
    assert all("third" in (a[2], b[2]) for a, b in asked)
    # A file that changed is compared with the others again.
    asked.clear()
    (tmp_path / "home" / "chunks" / "third.json").unlink()  # as the refusal says
    export(second, [*anions[:3], "CCCCC(=O)O", "[O-]C(=O)c1ccccc1"])
    assert _sdf_build(tmp_path / "home", first, second).exit_code == 0
    assert asked
    assert all({a[2], b[2]} == {"first", "second"} for a, b in asked)


def test_a_record_that_cannot_be_judged_does_not_hide_the_next_one(tmp_path):
    """Each ID was judged by its first record in a file. When that one had no
    structure, or both files carried "N/A" where the InChIKey goes, every record
    of the ID passed for the same compound, and aspirin came back as caffeine."""
    from rdkit import Chem

    def export(path, records):
        with Chem.SDWriter(str(path)) as writer:
            for number, smiles, key in records:
                mol = Chem.MolFromSmiles(smiles)
                mol.SetProp("_Name", str(number))
                if key:
                    mol.SetProp("PUBCHEM_IUPAC_INCHIKEY", key)
                writer.write(mol)
        return path

    empty_first = export(tmp_path / "a.sdf", [(1, "", None), (1, ASPIRIN, None)])
    other = export(tmp_path / "b.sdf", [(1, CAFFEINE, None)])
    result = _sdf_build(tmp_path / "home1", empty_first, other)
    assert result.exit_code == 1, result.output
    assert "ID 1 stands for different compounds" in result.output

    no_key = export(tmp_path / "c.sdf", [(7, "CCO", "N/A")])
    no_key_either = export(tmp_path / "d.sdf", [(7, "CCN", "N/A")])
    result = _sdf_build(tmp_path / "home2", no_key, no_key_either)
    assert result.exit_code == 1, result.output
    assert "ID 7 stands for different compounds: CCO in c and CCN in d" in result.output


def test_a_clash_among_chunks_that_record_no_file_is_not_a_crash(tmp_path):
    """Two chunks built by 2.2.0 clashing with each other, next to a new file with
    the same ID: IndexError instead of a refusal. And such a chunk is not one of
    PubChem's, whatever the message said about CIDs."""
    from pubchem_finder.index.builder import SdfSource, build_chunks
    from pubchem_finder.index.layout import IndexPaths
    from tests.conftest import write_sdf

    home = tmp_path / "home"
    paths = IndexPaths(home)
    for name, smiles in (("old_a", ASPIRIN), ("old_b", CAFFEINE)):
        path = write_sdf(
            tmp_path / f"{name}.sdf", [(1, smiles, "x", "title")], gz=False
        )
        build_chunks(paths, [SdfSource(name=name, local_path=path)])
        stats = paths.chunk_stats(name)
        data = json.loads(stats.read_text())
        del data["source"]  # as 2.2.0 wrote it
        stats.write_text(json.dumps(data))
    new = write_sdf(tmp_path / "new.sdf", [(1, PARACETAMOL, "x", "title")], gz=False)
    result = _sdf_build(home, new)
    assert result.exit_code == 1
    assert isinstance(result.exception, SystemExit)
    assert "ID 1 stands for different compounds" in result.output
    assert f"{PARACETAMOL} in new and {ASPIRIN} in old_a" in result.output
    assert "Give the records of new other IDs" in result.output
    assert "PubChem CID" not in result.output


def test_the_same_compound_with_and_without_stereo_is_still_merged(tmp_path):
    """Both records carry an InChIKey tag, as PubChem's do: only its first block
    (the skeleton) is compared."""
    from tests.conftest import write_sdf

    flat = "COc1ccc2cc(ccc2c1)C(C)C(=O)O"
    a = write_sdf(
        tmp_path / "a.sdf",
        [(156391, flat.replace("C(C)", "[C@H](C)"), "naproxen", "new")],
        gz=False,
    )
    b = write_sdf(tmp_path / "b.sdf", [(156391, flat, "flat", "new")], gz=False)
    result = _sdf_build(tmp_path / "home", a, b)
    assert result.exit_code == 0, result.output


def test_the_same_compound_in_two_files_is_still_merged(tmp_path, sample_sdfs):
    """An --sdf file that repeats compounds of another chunk, drawn another way:
    Kekule SMILES against PubChem's aromatic one, a charged form, no stereo."""
    from rdkit import Chem

    home = tmp_path / "home"
    mine = tmp_path / "mine.sdf"
    with Chem.SDWriter(str(mine)) as writer:
        for cid, smiles in (
            (2244, "CC(=O)OC1=CC=CC=C1C(O)=O"),  # aspirin, Kekule
            (338, "[O-]C(=O)c1ccccc1O"),  # salicylate, where the chunk has the acid
            (702, "OCC"),
        ):
            mol = Chem.MolFromSmiles(smiles)
            mol.SetProp("_Name", str(cid))
            mol.SetProp("PUBCHEM_SMILES", smiles)
            writer.write(mol)
    result = _sdf_build(home, sample_sdfs[0], mine)
    assert result.exit_code == 0, result.output
    assert "Index ready: 8 compounds from 2 chunk(s)" in result.output


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


@pytest.mark.parametrize(("limit", "truncated"), [(2, True), (3, False)])
def test_pubchem_exact_accepts_a_result_limit(pubchem_knowledge_api, limit, truncated):
    from tests.pubchem_fixtures import ASPIRIN, table

    pubchem_knowledge_api["identity"].respond(
        200, json=table(ASPIRIN, ASPIRIN | {"CID": 2245}, ASPIRIN | {"CID": 2246})
    )
    result = runner.invoke(
        app,
        ["pubchem", "exact", "CC(=O)Oc1ccccc1C(=O)O", "--max", str(limit), "--json"],
    )
    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["max_records"] == len(payload["hits"]) == limit
    assert payload["truncated"] is truncated
    assert (
        "MaxRecords"
        not in pubchem_knowledge_api["identity"].calls.last.request.url.params
    )


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
    no_names = runner.invoke(app, ["pubchem", "synonyms", "5", "--json"])
    assert no_names.exit_code == 0, no_names.output
    assert json.loads(no_names.output)["synonyms"] == []
    assert json.loads(no_names.output)["total"] == 0
    missing = runner.invoke(
        app, ["pubchem", "annotations", "2244", "--heading", "Toxicity"]
    )
    assert missing.exit_code == 1
    assert "no 'Toxicity' section" in missing.output
    bad = runner.invoke(app, ["pubchem", "similarity", "not smiles("])
    assert bad.exit_code == 1
    for args in (["synonyms", "8"], ["substructure", "not("], ["exact", "not("]):
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


def test_draw_reports_an_unwritable_output_path(tmp_path):
    missing = tmp_path / "missing" / "grid.png"
    result = runner.invoke(app, ["draw", "CCO", "-o", str(missing)])
    assert result.exit_code == 1
    assert result.output.startswith("error:")
    assert str(missing) in result.output


@pytest.mark.parametrize(
    ("command", "query"),
    [("search", "CCO"), ("substructure", "c1ccccc1")],
)
def test_search_reports_an_uncreatable_data_directory(
    tmp_path, monkeypatch, command, query
):
    blocker = tmp_path / "blocker"
    blocker.write_text("a file where the data directory should be")
    monkeypatch.setenv("PUBCHEM_FINDER_HOME", str(blocker / "home"))
    monkeypatch.setenv("PUBCHEM_FINDER_DEMO", "1")
    result = runner.invoke(app, [command, query])
    assert result.exit_code == 1
    assert "cannot create the data directory" in result.output
    assert "PUBCHEM_FINDER_HOME" in result.output


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
    assert (std[0]["removed"], std[0]["note"]) == ("[Na+]", "")
    assert (std[2]["removed"], std[2]["note"]) == ("", "")


def test_standardize_batch_says_what_was_left_out_and_why(tmp_path):
    src = tmp_path / "in.smi"
    src.write_text(
        "N.N.Cl[Pt]Cl cisplatin\nCC(C)(C)NCC(O)c1ccc(O)c(CO)c1.OS(=O)(=O)O.[Na+] two\n"
    )
    out = tmp_path / "standardize.csv"
    result = runner.invoke(app, ["standardize", "--input", str(src), "-o", str(out)])
    assert result.exit_code == 0, result.output
    std = list(csv.DictReader(out.open()))
    assert std[0]["removed"] == ""
    assert std[0]["note"].startswith("Contains Pt")
    assert std[1]["removed"] == "O=S(=O)(O)O.[Na+]"  # one SMILES, dot-separated


@pytest.mark.parametrize(
    "command", ["props", "standardize", "triage", "search", "substructure"]
)
@pytest.mark.parametrize(
    ("content", "expected"), [("nope( first\nnope( second\n", 1), ("", 0)]
)
def test_batch_exit_status_reports_no_successes(tmp_path, command, content, expected):
    src = tmp_path / "inputs.smi"
    src.write_text(content)
    output = tmp_path / "out.csv"
    result = runner.invoke(app, [command, "-i", str(src), "-o", str(output)])
    assert result.exit_code == expected, result.output
    rows = list(csv.DictReader(output.open()))
    assert len(rows) == (2 if content else 0)
    assert all(row["error"] for row in rows)
    assert f"{len(rows)} error(s)" in result.stderr


def test_batch_to_stdout_and_argument_checks(tmp_path):
    src = tmp_path / "in.csv"
    src.write_text("id,smiles\nethanol,CCO\n")
    result = runner.invoke(app, ["props", "-i", str(src)])
    assert result.exit_code == 0
    assert "ethanol,CCO,CCO," in result.output
    dst = tmp_path / "out.csv"
    result = runner.invoke(app, ["props", "-i", str(src), "-o", str(dst)])
    assert result.exit_code == 0
    header, record = dst.read_text().splitlines()
    columns = header.split(",")
    assert columns[3] == "molecular_formula"
    assert record.split(",")[3] == "C2H6O"
    both = runner.invoke(app, ["props", "CCO", "-i", str(src)])
    assert both.exit_code == 2
    neither = runner.invoke(app, ["props"])
    assert neither.exit_code == 2
    no_column = tmp_path / "bad.csv"
    no_column.write_text("name\nx\n")
    assert runner.invoke(app, ["props", "-i", str(no_column)]).exit_code == 1


@pytest.mark.parametrize("name", ["utf8_bom_smiles_first.csv", "utf8_bom_id_first.csv"])
@pytest.mark.parametrize("command", ["props", "standardize", "triage"])
def test_batch_utf8_bom_preserves_ids(tmp_path, name, command):
    src = Path(__file__).parents[1] / "fixtures" / "batch" / name
    out = tmp_path / "out.csv"
    result = runner.invoke(app, [command, "-i", str(src), "-o", str(out)])
    assert result.exit_code == 0, result.output
    with out.open(encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh))
    assert [(row["id"], row["smiles"], row["error"]) for row in rows] == [
        ("ethanol", "CCO", "")
    ]


@pytest.mark.parametrize("command", ["props", "standardize", "triage"])
def test_batch_stdin_utf8_bom_preserves_structure_and_id(command):
    result = runner.invoke(app, [command, "-i", "-"], input="\ufeffCCO ethanol\n")
    assert result.exit_code == 0, result.output
    rows = list(csv.DictReader(result.stdout.splitlines()))
    assert [(row["id"], row["smiles"], row["error"]) for row in rows] == [
        ("ethanol", "CCO", "")
    ]


@pytest.mark.parametrize(
    "command", ["props", "standardize", "triage", "search", "substructure"]
)
def test_batch_non_utf8_reports_error_without_emptying_output(tmp_path, command):
    src = Path(__file__).parents[1] / "fixtures" / "batch" / "windows_1252.csv"
    out = tmp_path / "out.csv"
    out.write_text("precious,data\n")
    result = runner.invoke(app, [command, "-i", str(src), "-o", str(out)])
    assert result.exit_code == 1
    assert "error:" in result.stderr
    assert "UTF-8" in result.stderr
    assert "save" in result.stderr.lower()
    assert str(src) in result.stderr
    assert "Traceback" not in result.output
    assert isinstance(result.exception, SystemExit)
    assert out.read_text() == "precious,data\n"


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


@pytest.mark.parametrize("bad", ["encoding", "field"])
def test_batch_late_input_error_preserves_previous_output(tmp_path, bad):
    import csv

    src = tmp_path / "late.csv"
    valid = b"id,smiles,comment\n" + b"r,CCO,ok\n" * 3000
    invalid = (
        b"caf\xe9,CCO,ok\n"
        if bad == "encoding"
        else b"r,CCO," + b"x" * (csv.field_size_limit() + 1) + b"\n"
    )
    src.write_bytes(valid + invalid)
    out = tmp_path / "keep.csv"
    out.write_bytes(b"previous results")
    result = runner.invoke(app, ["props", "-i", str(src), "-o", str(out)])
    assert result.exit_code == 1, result.output
    assert "error:" in result.stderr
    assert "Traceback" not in result.stderr
    assert out.read_bytes() == b"previous results"


def test_batch_reports_an_unwritable_output_before_the_run(tmp_path, monkeypatch):
    # A mistyped --output must not cost the whole batch: it is reported before
    # any row is processed, and names the file that was asked for.
    calls = []
    monkeypatch.setattr(service, "properties", calls.append)
    src = tmp_path / "in.smi"
    src.write_text("CCO\nCCN\n")
    missing = tmp_path / "missing-dir" / "out.csv"
    result = runner.invoke(app, ["props", "-i", str(src), "-o", str(missing)])
    assert result.exit_code == 1, result.output
    assert f"error: No such file or directory: {missing}" in result.stderr
    assert ".part" not in result.stderr
    assert calls == []


def _write_inputs(tmp_path):
    src = tmp_path / "in.smi"
    src.write_text("CCO first\nCCN second\nCCC third\n")
    out = tmp_path / "keep.csv"
    out.write_text("previous results\n")
    return src, out


def test_batch_keeps_the_previous_output_when_the_run_fails(tmp_path):
    """Issue #195: an existing output is only ever replaced by a complete result."""
    from pubchem_finder.cli import _run_batch

    src, out = _write_inputs(tmp_path)

    def breaks_on_the_second(row):
        if row.line == 2:
            raise RuntimeError("no space left")
        return [{"value": row.line}]

    with pytest.raises(RuntimeError):
        _run_batch(src, out, ["value"], breaks_on_the_second)
    assert out.read_text() == "previous results\n"
    assert sorted(p.name for p in tmp_path.iterdir()) == ["in.smi", "keep.csv"]


def _partial_files(folder):
    return sorted(p.name for p in folder.iterdir() if p.name.endswith(".part"))


def test_interrupted_batch_keeps_the_rows_done_and_the_previous_output(
    tmp_path, capsys
):
    from pubchem_finder.cli import _run_batch

    src, out = _write_inputs(tmp_path)

    def interrupted_on_the_third(row):
        if row.line == 3:
            raise KeyboardInterrupt
        return [{"value": row.line}]

    with pytest.raises(KeyboardInterrupt):
        _run_batch(src, out, ["value"], interrupted_on_the_third)
    (name,) = _partial_files(tmp_path)
    part = tmp_path / name
    assert name.startswith("keep.csv.")
    assert out.read_text() == "previous results\n"
    with part.open(newline="") as fh:
        assert [row["value"] for row in csv.DictReader(fh)] == ["1", "2"]
    assert f"the rows finished so far are in {part}" in capsys.readouterr().err


@pytest.mark.skipif(os.name == "nt", reason="POSIX links, file modes and pipes")
def test_batch_output_keeps_links_modes_and_special_files(tmp_path):
    import stat
    import threading

    src, _ = _write_inputs(tmp_path)
    real = tmp_path / "real.csv"
    real.write_text("old\n")
    real.chmod(0o640)
    link = tmp_path / "link.csv"
    link.symlink_to(real)
    result = runner.invoke(app, ["props", "-i", str(src), "-o", str(link)])
    assert result.exit_code == 0, result.output
    assert link.is_symlink()  # written through, not replaced by a file
    assert stat.S_IMODE(real.stat().st_mode) == 0o640
    assert len(real.read_text().splitlines()) == 4
    assert _partial_files(tmp_path) == []

    # A link to nothing yet: the file is created where the link points, as
    # open() does it.
    (tmp_path / "later.csv").symlink_to(tmp_path / "made.csv")
    result = runner.invoke(
        app, ["props", "-i", str(src), "-o", str(tmp_path / "later.csv")]
    )
    assert result.exit_code == 0, result.output
    assert (tmp_path / "later.csv").is_symlink()
    assert len((tmp_path / "made.csv").read_text().splitlines()) == 4

    # A pipe (or /dev/null) is written to directly: putting a regular file in its
    # place would break whatever reads it.
    fifo = tmp_path / "results.fifo"
    os.mkfifo(fifo)
    received: list[str] = []
    reader = threading.Thread(
        target=lambda: received.append(fifo.read_text()), daemon=True
    )
    reader.start()
    result = runner.invoke(app, ["props", "-i", str(src), "-o", str(fifo)])
    reader.join(timeout=30)
    assert result.exit_code == 0, result.output
    assert stat.S_ISFIFO(fifo.stat().st_mode)
    assert len(received[0].splitlines()) == 4


@pytest.mark.skipif(os.name == "nt", reason="names for open descriptors are POSIX")
@pytest.mark.parametrize("stdout_is", ["a pipe", "a file"])
def test_batch_writes_to_a_name_for_an_open_descriptor(tmp_path, stdout_is):
    # `-o /dev/stdout | head`, `-o >(gzip > out.csv.gz)`: the path names a pipe,
    # or a file the shell holds open. Looking for a file to replace behind it
    # found none ("No such file or directory: /dev/stdout") or, for a file that
    # was replaced once already, made one called "both.csv (deleted)".
    import subprocess
    import sys

    src, _ = _write_inputs(tmp_path)
    work = tmp_path / "work"
    work.mkdir()
    command = [sys.executable, "-m", "pubchem_finder", "props", "-i", str(src)]
    if stdout_is == "a pipe":
        result = subprocess.run(
            [*command, "-o", "/dev/stdout"],
            cwd=work,
            capture_output=True,
            text=True,
            timeout=120,
        )
        written = result.stdout
    else:
        with (work / "both.csv").open("w") as stdout:
            for _ in range(2):  # as in `{ props ...; props ...; } > both.csv`
                result = subprocess.run(
                    [*command, "-o", "/dev/stdout"],
                    cwd=work,
                    stdout=stdout,
                    stderr=subprocess.PIPE,
                    text=True,
                    timeout=120,
                )
                assert result.returncode == 0, result.stderr
        written = (work / "both.csv").read_text()
        assert [p.name for p in work.iterdir()] == ["both.csv"]
    assert result.returncode == 0, result.stderr
    assert "3 input(s), 3 row(s) written to /dev/stdout" in result.stderr
    ids = [line.split(",")[0] for line in written.splitlines()]
    # Opening /dev/stdout starts the file again on Linux and carries on where the
    # descriptor is on macOS, so two runs leave one copy of the rows or two.
    assert ids in (["id", "first", "second", "third"] * n for n in (1, 2))


@pytest.mark.parametrize(
    ("path", "descriptor"),
    [
        ("-", 0),
        ("/dev/stdin", 0),
        ("/dev/fd/0", 0),
        ("/dev/fd/7", 7),
        ("/proc/self/fd/12", 12),
        ("lib.smi", None),
        ("/dev/fd/x", None),
        ("/dev/stdin.smi", None),
    ],
)
def test_batch_input_names_for_open_descriptors(path, descriptor):
    from pubchem_finder.cli import _input_descriptor

    assert _input_descriptor(Path(path)) == descriptor


def test_batch_input_check_survives_a_descriptor_that_is_not_open(tmp_path):
    from pubchem_finder.cli import _same_file

    out = tmp_path / "out.csv"
    out.write_text("previous results\n")
    for descriptor in ("987654", "99999999999999999999"):  # the second overflows
        assert _same_file(Path(f"/dev/fd/{descriptor}"), out) is False
    assert _same_file(tmp_path / "missing.smi", out) is False
    assert _same_file(out, tmp_path / "missing.csv") is False


def test_batch_leaves_files_it_did_not_create_alone(tmp_path):
    # The partial file once had a fixed name, <output>.part, and whatever was
    # there was removed to make room: the rows an interrupted run had left, or
    # the input itself (`props -i lib.smi.part -o lib.smi` deleted its input,
    # and with it everything when the run then failed).
    src, _ = _write_inputs(tmp_path)
    content = src.read_text()
    for name in ("lib.smi.part", "out.csv.part"):
        (tmp_path / name).write_text(content)
    lib = tmp_path / "lib.smi"
    result = runner.invoke(
        app, ["props", "-i", str(tmp_path / "lib.smi.part"), "-o", str(lib)]
    )
    assert result.exit_code == 0, result.output
    assert len(lib.read_text().splitlines()) == 4
    result = runner.invoke(
        app, ["props", "-i", str(src), "-o", str(tmp_path / "out.csv")]
    )
    assert result.exit_code == 0, result.output
    assert len((tmp_path / "out.csv").read_text().splitlines()) == 4
    assert (tmp_path / "lib.smi.part").read_text() == content
    assert (tmp_path / "out.csv.part").read_text() == content
    assert _partial_files(tmp_path) == ["lib.smi.part", "out.csv.part"]


@pytest.mark.skipif(os.name == "nt", reason="POSIX links")
def test_batch_does_not_write_through_a_link_planted_as_its_partial_file(tmp_path):
    # In a directory others can write to, a partial file whose name can be
    # guessed could be put there first, as a link: the rows went to its target,
    # and the link was moved into place as the result.
    src, _ = _write_inputs(tmp_path)
    victim = tmp_path / "victim.txt"
    victim.write_text("not yours\n")
    out = tmp_path / "out.csv"
    (tmp_path / "out.csv.part").symlink_to(victim)
    result = runner.invoke(app, ["props", "-i", str(src), "-o", str(out)])
    assert result.exit_code == 0, result.output
    assert victim.read_text() == "not yours\n"
    assert out.is_file()
    assert not out.is_symlink()
    assert len(out.read_text().splitlines()) == 4
    assert (tmp_path / "out.csv.part").is_symlink()


@pytest.mark.skipif(os.name == "nt", reason="POSIX file modes")
def test_batch_keeps_the_rows_of_a_private_output_private(tmp_path):
    import stat

    from pubchem_finder.cli import _run_batch

    src, out = _write_inputs(tmp_path)
    out.chmod(0o600)
    seen = []

    def look(row):
        seen.extend(
            stat.S_IMODE((tmp_path / name).stat().st_mode)
            for name in _partial_files(tmp_path)
        )
        return [{"value": row.line}]

    _run_batch(src, out, ["value"], look)
    assert seen == [0o600, 0o600, 0o600]
    assert stat.S_IMODE(out.stat().st_mode) == 0o600
    assert len(out.read_text().splitlines()) == 4


def _fails_on_the_second(row):
    if row.line == 2:
        raise RuntimeError("no space left")
    return [{"value": row.line}]


@pytest.mark.parametrize("obstacle", ["no room beside it", "cannot be replaced"])
def test_batch_writes_in_place_an_output_it_cannot_replace(
    tmp_path, monkeypatch, obstacle
):
    # Replacing a file needs more than writing it does: a place for the partial
    # file beside it (not in a read-only directory) and leave to rename over it
    # (not for a file mounted into a container, nor in a directory that lets
    # only a file's owner rename it). Such an output was refused, or the run
    # failed after its last row. Its rows are staged all the same and copied
    # into it at the end.
    import errno
    import tempfile

    from pubchem_finder import cli

    src, out = _write_inputs(tmp_path)
    before = out.stat()
    staging = tmp_path / "system-tmp"
    staging.mkdir()
    monkeypatch.setattr(tempfile, "tempdir", str(staging))
    if obstacle == "no room beside it":

        def refuse(target, previous):
            raise PermissionError(errno.EACCES, "Permission denied", str(target))

        monkeypatch.setattr(cli, "_stage_beside", refuse)
    else:

        def busy(src, dst):
            raise OSError(errno.EBUSY, "Device or resource busy", str(dst))

        monkeypatch.setattr(cli.os, "replace", busy)

    with pytest.raises(RuntimeError):
        cli._run_batch(src, out, ["value"], _fails_on_the_second)
    assert out.read_text() == "previous results\n"

    result = runner.invoke(app, ["props", "-i", str(src), "-o", str(out)])
    assert result.exit_code == 0, result.output
    assert len(out.read_text().splitlines()) == 4
    after = out.stat()
    assert (after.st_dev, after.st_ino) == (before.st_dev, before.st_ino)
    assert _partial_files(tmp_path) == []
    assert list(staging.iterdir()) == []


def test_batch_says_where_the_rows_are_when_the_output_cannot_be_written(
    tmp_path, monkeypatch
):
    import errno

    from pubchem_finder import cli

    src, _ = _write_inputs(tmp_path)
    out = tmp_path / "new.csv"

    def busy(src, dst):
        raise OSError(errno.EBUSY, "Device or resource busy", str(dst))

    monkeypatch.setattr(cli.os, "replace", busy)
    result = runner.invoke(app, ["props", "-i", str(src), "-o", str(out)])
    assert result.exit_code == 1, result.output
    (name,) = _partial_files(tmp_path)
    assert (
        f"error: could not write {out} (Device or resource busy); "
        f"the results are in {tmp_path / name}"
    ) in result.stderr.replace("\n", "")
    assert len((tmp_path / name).read_text().splitlines()) == 4


@pytest.mark.skipif(os.name == "nt", reason="POSIX hard links")
def test_batch_keeps_the_other_names_of_its_output(tmp_path):
    # Renaming a new file over one name of a file leaves its other names with
    # the old content.
    src, out = _write_inputs(tmp_path)
    other = tmp_path / "other-name.csv"
    os.link(out, other)
    result = runner.invoke(app, ["props", "-i", str(src), "-o", str(out)])
    assert result.exit_code == 0, result.output
    assert len(other.read_text().splitlines()) == 4
    assert out.stat().st_ino == other.stat().st_ino
    assert _partial_files(tmp_path) == []


def _file_permissions_apply(tmp_path):
    """False for root, which may read and write whatever the modes say."""
    probe = tmp_path / "permission-probe"
    probe.mkdir()
    probe.chmod(0)
    try:
        return not os.access(probe, os.R_OK | os.W_OK | os.X_OK)
    finally:
        probe.chmod(0o755)
        probe.rmdir()


@pytest.mark.skipif(os.name == "nt", reason="POSIX directory modes")
def test_batch_writes_a_file_in_a_directory_it_cannot_write(tmp_path):
    if not _file_permissions_apply(tmp_path):
        pytest.skip("this user is not held to file modes")
    src, _ = _write_inputs(tmp_path)
    results = tmp_path / "results"
    results.mkdir()
    out = results / "out.csv"
    out.write_text("previous results\n")
    results.chmod(0o555)
    try:
        result = runner.invoke(app, ["props", "-i", str(src), "-o", str(out)])
        assert result.exit_code == 0, result.output
        assert len(out.read_text().splitlines()) == 4
        # A file that does not exist yet cannot be made there, which is said
        # before the first row is processed.
        missing = results / "new.csv"
        result = runner.invoke(app, ["props", "-i", str(src), "-o", str(missing)])
        assert result.exit_code == 1, result.output
        assert f"error: Permission denied: {missing}" in result.stderr
    finally:
        results.chmod(0o755)


@pytest.mark.skipif(os.name == "nt", reason="POSIX links and directory modes")
def test_batch_does_not_follow_a_link_the_system_refuses(tmp_path):
    # A link the kernel will not follow (one planted in /tmp by another user,
    # where fs.protected_symlinks is on; here, one through a directory that may
    # not be entered) was resolved by reading it, and its target was replaced
    # with the results.
    if not _file_permissions_apply(tmp_path):
        pytest.skip("this user is not held to file modes")
    src, _ = _write_inputs(tmp_path)
    shared = tmp_path / "shared"
    (shared / "closed").mkdir(parents=True)
    victim = shared / "victim.txt"
    victim.write_text("not yours\n")
    out = shared / "out.csv"
    out.symlink_to("closed/../victim.txt")
    (shared / "closed").chmod(0)
    try:
        for source in (["-i", str(src)], ["-i", "-"]):
            result = runner.invoke(
                app, ["props", *source, "-o", str(out)], input="CCO first\n"
            )
            assert result.exit_code == 1, result.output
            assert f"error: Permission denied: {out}" in result.stderr
            assert "Traceback" not in result.output
            assert isinstance(result.exception, SystemExit)
    finally:
        (shared / "closed").chmod(0o755)
    assert victim.read_text() == "not yours\n"
    assert out.is_symlink()


def test_batch_asks_the_system_before_following_a_link(tmp_path, monkeypatch):
    # The same refusal, for users who are not held to file modes (root in a
    # container): what the system says about the link decides.
    from pubchem_finder import cli

    real = tmp_path / "real.csv"
    real.write_text("previous results\n")
    link = tmp_path / "link.csv"
    try:
        link.symlink_to(real)
    except OSError:
        pytest.skip("no permission to create links here")
    plan = cli._replaceable(link)
    assert plan is not None
    assert plan[1] is not None
    assert plan[0] == Path(os.path.realpath(real))  # the file, not the link
    assert plan[1].st_ino == real.stat().st_ino
    plain = cli._replaceable(real)
    assert plain is not None
    assert plain[0] == real
    assert cli._replaceable(tmp_path / "new.csv") == (tmp_path / "new.csv", None)
    assert cli._replaceable(tmp_path) is None  # not a file: open() reports it

    def refused(path, *args, **kwargs):
        raise PermissionError(13, "Permission denied", str(path))

    monkeypatch.setattr(cli.os, "stat", refused)
    assert cli._replaceable(link) is None


def _library(tmp_path, n=1500):
    """A .smi file larger than a pipe's buffer (64 KiB on Linux)."""
    lib = tmp_path / "lib.smi"
    lib.write_text(
        "".join(f"C{'C' * (i % 5)}O compound-{i:05d}-{'x' * 80}\n" for i in range(n))
    )
    return lib


@pytest.mark.parametrize(
    "spelling",
    [
        "-",
        pytest.param(
            "/dev/stdin",
            marks=pytest.mark.skipif(os.name == "nt", reason="no /dev/stdin"),
        ),
    ],
)
def test_batch_refuses_to_write_over_redirected_standard_input(tmp_path, spelling):
    """Issue #195: `props -i - -o lib.smi < lib.smi` left 179 of 5,000 rows, exit 0."""
    import subprocess
    import sys

    lib = _library(tmp_path)
    before = lib.read_bytes()
    command = [sys.executable, "-m", "pubchem_finder", "props"]
    with lib.open("rb") as stdin:
        result = subprocess.run(
            [*command, "-i", spelling, "-o", str(lib)],
            stdin=stdin,
            capture_output=True,
            text=True,
            timeout=120,
        )
    assert result.returncode == 2, result.stderr
    assert "is the input file" in result.stderr
    assert lib.read_bytes() == before


@pytest.mark.skipif(
    os.name == "nt", reason="Windows cannot replace a file another process reads"
)
def test_batch_reads_piped_input_to_its_end_before_replacing_its_source(tmp_path):
    """`cat lib.smi | props -i - -o lib.smi` cannot be recognised as writing over
    the input, so the results replace the file only once every row was read:
    emptying it first lost 1,522 of 5,000 rows, with exit 0."""
    import subprocess
    import sys

    lib = _library(tmp_path)
    feed = "import shutil, sys; shutil.copyfileobj(open(sys.argv[1], 'rb'), sys.stdout.buffer)"
    feeder = subprocess.Popen(
        [sys.executable, "-c", feed, str(lib)], stdout=subprocess.PIPE
    )
    result = subprocess.run(
        [sys.executable, "-m", "pubchem_finder", "props", "-i", "-", "-o", str(lib)],
        stdin=feeder.stdout,
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert feeder.stdout is not None
    feeder.stdout.close()
    assert feeder.wait(timeout=30) == 0
    assert result.returncode == 0, result.stderr
    assert "1500 input(s), 1500 row(s) written" in result.stderr
    with lib.open(newline="") as fh:
        rows = list(csv.DictReader(fh))
    assert len(rows) == 1500
    assert not any(row["error"] for row in rows)


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
    assert "note:" not in result.stderr
    rows = list(csv.DictReader(out.open()))
    aspirin = [r for r in rows if r["id"] == "aspirin"]
    assert aspirin[0]["rank"] == "1"
    assert aspirin[0]["cid"] == "2244"
    assert any(r["id"] == "ethanol" for r in rows)

    # Issue #191 for a file of queries: the rows of an input are its best --top
    # hits, and nothing said so when more reach the threshold. Six compounds
    # reach 0.2 for aspirin here and one for ethanol.
    def cut_at(top):
        args = ["search", "-i", str(src), "-o", str(out), "--home", str(built_home)]
        return runner.invoke(app, [*args, "--threshold", "0.2", "--top", top])

    cut = cut_at("2")
    assert cut.exit_code == 0, cut.output
    assert (
        "note: 1 input(s) have more hits than --top 2 (aspirin): their rows are "
        "the best 2, not all; raise --top (up to 1000) or --threshold"
    ) in cut.stderr
    assert len([r for r in csv.DictReader(out.open()) if r["id"] == "aspirin"]) == 2
    longest = cut_at("1000")
    assert longest.exit_code == 0, longest.output
    assert "note:" not in longest.stderr
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


def test_resolve_home_treats_blank_as_default(monkeypatch):
    from pubchem_finder.config import default_home, resolve_home

    monkeypatch.delenv("PUBCHEM_FINDER_HOME", raising=False)
    assert resolve_home("") == default_home()
    assert resolve_home("   ") == default_home()


# --- `--home` on the live PubChem commands moves their response cache (#161) --


def test_lookup_home_puts_its_cache_there(tmp_path, pubchem_api, monkeypatch):
    """A lookup run with --home caches its response under that data directory."""
    monkeypatch.setattr(service, "_clients", {})
    monkeypatch.setenv("PUBCHEM_FINDER_CACHE", "1")
    monkeypatch.setattr("pubchem_finder.data.pug_rest.time.sleep", lambda s: None)
    home = tmp_path / "home"
    result = runner.invoke(app, ["lookup", "caffeine", "--home", str(home)])
    assert result.exit_code == 0, result.output
    assert "CID 2519  Caffeine" in result.output
    assert (home / "cache" / "pubchem-http.sqlite").exists()


def test_pubchem_commands_share_the_home_cache_they_report(
    tmp_path, pubchem_knowledge_api, monkeypatch
):
    """`cache info --home DIR` sees what `pubchem ... --home DIR` cached."""
    monkeypatch.setattr(service, "_clients", {})
    monkeypatch.setenv("PUBCHEM_FINDER_CACHE", "1")
    monkeypatch.setattr("pubchem_finder.data.pug_rest.time.sleep", lambda s: None)
    home = tmp_path / "home"
    synonyms = runner.invoke(app, ["pubchem", "synonyms", "2244", "--home", str(home)])
    assert synonyms.exit_code == 0, synonyms.output
    described = runner.invoke(app, ["pubchem", "describe", "2244", "--home", str(home)])
    assert described.exit_code == 0, described.output
    info = runner.invoke(app, ["cache", "info", "--home", str(home)])
    assert info.exit_code == 0, info.output
    assert str(home.resolve() / "cache" / "pubchem-http.sqlite") in info.output
    assert "0 responses" not in info.output


def test_empty_home_option_is_rejected(monkeypatch, tmp_path):
    monkeypatch.chdir(tmp_path)
    monkeypatch.delenv("PUBCHEM_FINDER_HOME", raising=False)
    result = runner.invoke(app, ["search", "CCO", "--home", ""])
    assert result.exit_code == 2
    assert not (tmp_path / "demo").exists()


@pytest.mark.parametrize(
    "mode", ["--demo", "--with-substructure", "--without-substructure"]
)
@pytest.mark.parametrize(("option", "value"), [("--fp-size", "512"), ("--radius", "3")])
def test_fingerprint_options_are_refused_without_chunk_sources(
    tmp_path, mode, option, value
):
    home = tmp_path / "home"
    result = runner.invoke(
        app, ["index", "build", mode, option, value, "--home", str(home)]
    )
    assert result.exit_code == 2, result.output
    assert "applies only when chunks are built" in result.output
    assert not home.exists()


@pytest.mark.parametrize("finalized", [False, True])
def test_extending_an_index_keeps_its_fingerprint_settings(
    tmp_path, sample_sdfs, finalized
):
    from pubchem_finder.index.layout import IndexPaths

    home = tmp_path / "home"
    common = ["index", "build", "--home", str(home), "--workers", "1"]
    first = runner.invoke(
        app,
        [*common, "--sdf", str(sample_sdfs[0]), "--radius", "3", "--fp-size", "512"],
    )
    assert first.exit_code == 0, first.output
    paths = IndexPaths(home)
    before = paths.load_manifest()
    assert before.fp_params == {"radius": 3, "fpSize": 512}
    if not finalized:
        paths.manifest.unlink()

    extended = runner.invoke(app, [*common, "--sdf", str(sample_sdfs[1])])
    assert extended.exit_code == 0, extended.output
    after = paths.load_manifest()
    assert after.fp_params == before.fp_params
    assert after.n_compounds > before.n_compounds

    matching = runner.invoke(
        app, [*common, "--sdf", str(sample_sdfs[1]), "--fp-size", "512"]
    )
    assert matching.exit_code == 0, matching.output
    mismatched = runner.invoke(
        app, [*common, "--sdf", str(sample_sdfs[1]), "--radius", "2"]
    )
    assert mismatched.exit_code == 1, mismatched.output
    assert "Existing chunks" in mismatched.output
    assert paths.load_manifest() == after
