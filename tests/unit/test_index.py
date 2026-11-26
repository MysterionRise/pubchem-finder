"""Tests for building and searching the local index."""

import json
import os
from pathlib import Path

import numpy as np
import pytest
import tables as tb

from pubchem_finder import hints, service
from pubchem_finder.core.exceptions import (
    IndexBuildError,
    IndexIncompatible,
    IndexNotBuilt,
    InvalidInput,
    MoleculeValidationError,
    PubChemFinderError,
)
from pubchem_finder.index import builder as builder_mod
from pubchem_finder.index import engine as engine_mod
from pubchem_finder.index.builder import (
    SdfSource,
    build_chunks,
    finalize,
    h5_versions,
    runtime_versions,
)
from pubchem_finder.index.engine import LocalIndex, choose_in_memory
from pubchem_finder.index.layout import IndexPaths, chunk_name
from tests.conftest import (
    CHUNK_A,
    CHUNK_NAMES,
    IMINES,
    N_GOOD,
    pubchem_block,
    write_sdf,
)

ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"


def sources(paths):
    return [SdfSource(name=chunk_name(p), local_path=p) for p in paths]


@pytest.mark.parametrize("name", ["Library.SDF", "Library.Sd.GZ", "Library.MOL.gz"])
def test_chunk_name_strips_case_insensitive_suffixes(name):
    assert chunk_name(name) == "Library"


def test_uppercase_suffix_reuses_chunk_from_previous_release(tmp_path, sample_sdfs):
    sdf = tmp_path / "Library.SDF.GZ"
    sdf.write_bytes(sample_sdfs[0].read_bytes())
    paths = IndexPaths(tmp_path / "home")
    build_chunks(paths, [SdfSource(name=sdf.name, local_path=sdf)], workers=1)
    before = finalize(paths)

    build_chunks(paths, sources([sdf]), workers=1)
    after = finalize(paths)

    assert after.chunks == before.chunks == [sdf.name]
    assert after.n_compounds == before.n_compounds > 0


def test_legacy_chunk_is_reused_when_file_is_typed_in_another_case(
    tmp_path, sample_sdfs
):
    sdf = tmp_path / "Library.SDF.GZ"
    sdf.write_bytes(sample_sdfs[0].read_bytes())
    alias = tmp_path / "library.SDF.GZ"
    if not alias.exists():
        os.link(sdf, alias)
    paths = IndexPaths(tmp_path / "home")
    build_chunks(paths, [SdfSource(name=sdf.name, local_path=sdf)], workers=1)
    before = finalize(paths)

    build_chunks(paths, sources([alias]), workers=1)
    after = finalize(paths)

    assert after.chunks == before.chunks == [sdf.name]
    assert after.n_compounds == before.n_compounds


def test_legacy_chunk_is_reused_when_new_name_belongs_to_another_file(
    tmp_path, sample_sdfs
):
    modern = tmp_path / "a" / "LIB.sdf.gz"
    legacy = tmp_path / "b" / "LIB.SDF.GZ"
    for path, source in zip((modern, legacy), sample_sdfs, strict=False):
        path.parent.mkdir()
        path.write_bytes(source.read_bytes())
    paths = IndexPaths(tmp_path / "home")
    build_chunks(
        paths,
        [*sources([modern]), SdfSource(name=legacy.name, local_path=legacy)],
        workers=1,
    )
    before = finalize(paths)

    report = build_chunks(paths, sources([legacy]), workers=1)
    after = finalize(paths)

    assert report.skipped == [legacy.name]
    assert after.chunks == before.chunks
    assert after.n_compounds == before.n_compounds


def test_finalized_layout_and_manifest(built_home):
    paths = IndexPaths(built_home)
    manifest = paths.load_manifest()
    assert manifest is not None
    assert manifest.n_compounds == N_GOOD
    assert manifest.n_failed == 3  # two in chunk A, one in chunk C
    assert manifest.chunks == list(CHUNK_NAMES[:2])  # the all-invalid chunk is left out
    assert manifest.fp_type == "Morgan"
    assert manifest.fp_params == {"radius": 2, "fpSize": 1024}
    assert manifest.fingerprint == "Morgan fpSize=1024 radius=2"
    assert manifest.memory_bytes() == N_GOOD * (1024 // 64 + 2) * 8
    for name in CHUNK_NAMES:
        assert paths.chunk_stats(name).exists()
        assert paths.chunk_parquet(name).exists()
    assert not paths.chunk_h5(CHUNK_NAMES[2]).exists()
    assert not list(paths.chunks.glob("*.tmp"))
    assert json.loads(paths.manifest.read_text())["n_compounds"] == N_GOOD


def test_rebuild_skips_completed_chunks(built_home, sample_sdfs):
    report = build_chunks(IndexPaths(built_home), sources(sample_sdfs), workers=1)
    assert report.built == []
    assert report.skipped == list(CHUNK_NAMES)


def test_mismatched_fingerprint_params_rejected(built_home, sample_sdfs):
    with pytest.raises(IndexBuildError, match="Existing chunks"):
        build_chunks(
            IndexPaths(built_home),
            sources(sample_sdfs),
            fp_params={"radius": 3, "fpSize": 2048},
        )


def test_single_chunk_finalize_and_downloaded_file_cleanup(tmp_path, sample_sdfs):
    copy = tmp_path / sample_sdfs[1].name
    copy.write_bytes(sample_sdfs[1].read_bytes())
    paths = IndexPaths(tmp_path / "home")
    src = SdfSource(name=chunk_name(copy), fetch=lambda: copy, delete_after=True)
    events: list[str] = []
    build_chunks(paths, [src], workers=1, on_event=events.append)
    assert not copy.exists()  # downloaded SDF removed once its chunk is complete
    manifest = finalize(paths, on_event=events.append)
    assert manifest.n_compounds == 6
    assert any("built" in e for e in events)
    assert any("index ready" in e for e in events)


def test_finalize_without_chunks(tmp_path):
    with pytest.raises(IndexBuildError, match="nothing to finalize"):
        finalize(IndexPaths(tmp_path))


def test_source_without_path_or_fetch():
    with pytest.raises(IndexBuildError):
        SdfSource(name="x").path()


def test_parallel_build_matches_serial(tmp_path, sample_sdfs):
    paths = IndexPaths(tmp_path)
    report = build_chunks(paths, sources(sample_sdfs[:2]), workers=2)
    assert sorted(s.name for s in report.built) == list(CHUNK_NAMES[:2])
    assert finalize(paths).n_compounds == N_GOOD


def test_search_finds_aspirin_first(built_home):
    result = service.similarity_search(ASPIRIN, threshold=0.3, home=built_home)
    assert result.n_indexed == N_GOOD
    assert result.fingerprint.startswith("Morgan")
    top = result.hits[0]
    assert (top.cid, top.similarity) == (2244, 1.0)
    assert top.iupac_name == "2-acetyloxybenzoic acid"
    assert top.molecular_formula == "C9H8O4"
    assert top.pubchem_url.endswith("/2244")
    cids = [h.cid for h in result.hits]
    assert 338 in cids  # salicylic acid
    assert 2519 not in cids  # caffeine is not similar
    sims = [h.similarity for h in result.hits]
    assert sims == sorted(sims, reverse=True)


def test_search_merged_says_whether_more_reach_the_threshold(built_home):
    from rdkit import Chem

    index = LocalIndex(built_home)
    query = Chem.MolFromSmiles(ASPIRIN)
    everything = index.search_merged(query, 0.1, 1000)
    n = len(everything.pairs)
    assert n > 2
    assert everything.truncated is False
    assert index.search_merged(query, 0.1, n).truncated is False
    cut = index.search_merged(query, 0.1, n - 1)
    assert cut.truncated is True
    assert cut.pairs == everything.pairs[: n - 1]
    assert index.search_merged(query, 0.1, 1).truncated is True
    assert index.search_merged(query, 0.99, 1).truncated is False  # aspirin alone


@pytest.mark.parametrize("tag", [True, False])
def test_a_compound_scores_one_against_the_smiles_stored_for_it(tmp_path, tag):
    """Issue #187: 0.8% of PubChem (compounds with an imine or amidine N-H) were
    fingerprinted with an extra hydrogen atom and did not match themselves:
    metformin 0.68, benzamidine 0.67, guanidine 0.5."""
    library = tmp_path / "imines.sdf"
    library.write_text(
        "".join(pubchem_block(cid, smiles, tag) for cid, smiles in IMINES)
    )
    paths = IndexPaths(tmp_path / "home")
    report = build_chunks(paths, sources([library]), workers=1)
    finalize(paths)
    assert report.built[0].fingerprint_of == "smiles"
    stored = service.get_index(paths.home).metadata([cid for cid, _ in IMINES])
    for cid, smiles in IMINES:
        for query in (stored[cid]["smiles"], smiles):
            found = service.similarity_search(query, threshold=1.0, home=paths.home)
            assert [(hit.cid, hit.similarity) for hit in found.hits] == [(cid, 1.0)]


def test_every_indexed_compound_finds_itself(built_home):
    from tests.conftest import CHUNK_B

    cids = [cid for cid, *_ in CHUNK_A + CHUNK_B]
    stored = service.get_index(built_home).metadata(cids)
    assert len(stored) == N_GOOD
    for cid, row in stored.items():
        found = service.similarity_search(row["smiles"], threshold=1.0, home=built_home)
        assert cid in [hit.cid for hit in found.hits], row["smiles"]


def test_threshold_and_top_k(built_home):
    strict = service.similarity_search(ASPIRIN, threshold=0.99, home=built_home)
    assert [h.cid for h in strict.hits] == [2244]
    one = service.similarity_search(ASPIRIN, threshold=0.1, top_k=1, home=built_home)
    assert len(one.hits) == 1


@pytest.mark.parametrize("metric", ["dice", "cosine"])
def test_other_metrics(built_home, metric):
    result = service.similarity_search(
        ASPIRIN, threshold=0.3, metric=metric, home=built_home
    )
    assert result.hits[0].cid == 2244
    assert result.metric == metric


def test_disk_mode_matches_memory_mode(built_home):
    mem = LocalIndex(built_home, mode="memory")
    disk = LocalIndex(built_home, mode="disk")
    from rdkit import Chem

    query = Chem.MolFromSmiles(ASPIRIN)
    assert mem.search(query, 0.2, 10) == disk.search(query, 0.2, 10)
    assert mem.in_memory is True
    assert disk.in_memory is False


def test_search_input_validation(built_home):
    with pytest.raises(MoleculeValidationError):
        service.similarity_search("not smiles(", home=built_home)
    with pytest.raises(PubChemFinderError, match="threshold"):
        service.similarity_search(ASPIRIN, threshold=0, home=built_home)
    with pytest.raises(PubChemFinderError, match="top_k"):
        service.similarity_search(ASPIRIN, top_k=0, home=built_home)
    with pytest.raises(PubChemFinderError, match="metric"):
        service.similarity_search(ASPIRIN, metric="euclid", home=built_home)  # type: ignore[arg-type]


def test_search_without_index(tmp_path):
    with pytest.raises(IndexNotBuilt, match="index build"):
        service.similarity_search(ASPIRIN, home=tmp_path)


def test_metadata_of_empty_list(built_home):
    assert LocalIndex(built_home).metadata([]) == {}


def test_status_not_built_and_ready(tmp_path, built_home):
    missing = service.index_status(tmp_path)
    assert missing.state == "not_built"
    assert "index build" in (missing.hint or "")
    ready = service.index_status(built_home)
    assert ready.state == "ready"
    assert ready.n_compounds == N_GOOD
    assert ready.chunks_in_index == 2
    assert ready.chunks_built == 3
    assert ready.size_bytes > 0
    assert ready.built_at


def test_status_offers_to_finish_an_interrupted_build(tmp_path, built_home):
    first = sorted(IndexPaths(built_home).chunks.glob("*.json"))[0]
    IndexPaths(tmp_path).chunks.mkdir(parents=True)
    (IndexPaths(tmp_path).chunks / first.name).write_bytes(first.read_bytes())
    status = service.index_status(tmp_path)
    assert status.state == "not_built"
    assert status.chunks_built == 1
    assert "stopped after 1 chunk" in (status.hint or "")
    assert "323 MB" not in (status.hint or "")


def test_status_reports_loaded_mode(built_home):
    service.similarity_search(ASPIRIN, home=built_home)
    status = service.index_status(built_home)
    assert status.loaded is True
    assert status.mode in ("memory", "disk")


def test_choose_in_memory(monkeypatch, built_home):
    manifest = IndexPaths(built_home).load_manifest()
    assert manifest is not None
    assert choose_in_memory(manifest.memory_bytes(), "memory") is True
    assert choose_in_memory(manifest.memory_bytes(), "disk") is False
    monkeypatch.setattr(engine_mod, "available_memory", lambda: 10)
    assert choose_in_memory(manifest.memory_bytes(), "auto") is False
    monkeypatch.setattr(engine_mod, "available_memory", lambda: None)
    assert choose_in_memory(manifest.memory_bytes(), "auto") is True


def test_mode_from_environment(monkeypatch, built_home):
    monkeypatch.setenv("PUBCHEM_FINDER_INDEX_MODE", "disk")
    assert LocalIndex(built_home).mode == "disk"
    monkeypatch.setenv("PUBCHEM_FINDER_INDEX_MODE", "auto")
    assert LocalIndex(built_home).mode == "auto"
    monkeypatch.setenv("PUBCHEM_FINDER_INDEX_MODE", "bogus")
    assert LocalIndex(built_home, mode="memory").mode == "memory"
    with pytest.raises(
        InvalidInput, match="PUBCHEM_FINDER_INDEX_MODE must be auto, memory, or disk"
    ):
        LocalIndex(built_home)


def test_available_memory_is_positive_or_unknown():
    value = engine_mod.available_memory()
    assert value is None or value > 0


def test_finalize_skips_when_up_to_date_unless_forced(tmp_path, sample_sdfs):
    paths = IndexPaths(tmp_path)
    build_chunks(paths, sources(sample_sdfs[:1]), workers=1)
    first = finalize(paths)
    events: list[str] = []
    again = finalize(paths, on_event=events.append)
    assert again.built_at == first.built_at
    assert events == [f"index up to date: {first.n_compounds:,} compounds"]
    forced = finalize(paths, force=True)
    assert forced.built_at > first.built_at
    # adding a chunk triggers a real merge
    build_chunks(paths, sources(sample_sdfs[1:2]), workers=1)
    assert finalize(paths).n_compounds == N_GOOD


def test_similarities_are_clean_four_decimal_floats(built_home):
    result = service.similarity_search(ASPIRIN, threshold=0.3, home=built_home)
    for hit in result.hits:
        assert hit.similarity == round(hit.similarity, 4)
        assert len(repr(hit.similarity)) <= 6


@pytest.mark.parametrize("pattern", [False, True])
def test_sort_matches_fpsim2(tmp_path, sample_sdfs, pattern):
    """Our popcount sort writes the same file FPSim2's sort_db_file does."""
    import shutil

    import tables as tb
    from FPSim2.io.backends.pytables import sort_db_file

    from pubchem_finder.index.builder import build_pattern_chunk, sort_fpsim2_file

    paths = IndexPaths(tmp_path / "home")
    build_chunks(paths, sources(sample_sdfs[1:2]), workers=1)
    name = CHUNK_NAMES[1]
    src = paths.chunk_h5(name)
    if pattern:
        build_pattern_chunk(str(paths.home), name)
        src = paths.pattern_h5(name)
    ours, ref = tmp_path / "ours.h5", tmp_path / "ref.h5"
    sort_fpsim2_file(src, ours)
    shutil.copyfile(src, ref)
    sort_db_file(str(ref))
    with tb.open_file(str(ours)) as a, tb.open_file(str(ref)) as b:
        assert [a.root.config[i] for i in range(4)] == [
            b.root.config[i] for i in range(4)
        ]
        assert a.root.config[4] == [
            (int(v), (int(lo), int(hi))) for v, (lo, hi) in b.root.config[4]
        ]
        assert (a.root.fps[:] == b.root.fps[:]).all()


def test_fpsim2_version_warning_stays_off_stdout(built_home, monkeypatch, capsys):
    """FPSim2 print()s when RDKit versions differ; stdout is the MCP transport."""
    import rdkit
    from rdkit import Chem

    monkeypatch.setattr(rdkit, "__version__", "2000.01.1")
    LocalIndex(built_home, mode="memory").search(Chem.MolFromSmiles(ASPIRIN), 0.5, 5)
    captured = capsys.readouterr()
    assert captured.out == ""
    assert "RDKit version" in captured.err


def test_identical_fingerprints_break_tie_by_cid(tmp_path):
    """Two compounds with identical fingerprints tie on score: top_k must cut
    after IDs break the tie, so the lower CID wins whoever FPSim2 returned."""
    from rdkit import Chem

    aspirin = "CC(=O)Oc1ccccc1C(=O)O"
    # The 11C twin of aspirin has the same Morgan fingerprint; store it first so
    # FPSim2's own order favours it when the search cuts at top_k.
    sdf = write_sdf(
        tmp_path / "Compound_000000001_000500000.sdf.gz",
        [
            (450661, aspirin, "11C-labelled aspirin", "new"),
            (2244, aspirin, "2-acetyloxybenzoic acid", "new"),
        ],
    )
    paths = IndexPaths(tmp_path)
    build_chunks(paths, sources([sdf]), workers=1)
    finalize(paths)
    index = LocalIndex(tmp_path, mode="memory")
    query = Chem.MolFromSmiles(aspirin)
    assert index.search(query, 0.5, 1) == [(2244, 1.0)]
    assert index.search(query, 0.5, 2) == [(2244, 1.0), (450661, 1.0)]


def test_search_breaks_ties_by_cid(built_home, monkeypatch):
    """FPSim2 returns equal scores in thread order; ties must come lowest CID first."""
    import numpy as np
    from rdkit import Chem

    hits = np.array(
        [(450661, 1.0), (2244, 1.0), (4133, 0.7), (338, 0.7), (3672, 0.9)],
        dtype=[("mol_id", "<u4"), ("coeff", "<f4")],
    )

    class StubEngine:
        def top_k(self, *args, **kwargs):
            return hits

    index = LocalIndex(built_home, mode="memory")
    monkeypatch.setattr(index, "_ensure_loaded", lambda: (StubEngine(), index.merged))
    index._in_memory = True
    assert index.search(Chem.MolFromSmiles(ASPIRIN), 0.5, 5) == [
        (2244, 1.0),
        (450661, 1.0),
        (3672, 0.9),
        (338, 0.7),
        (4133, 0.7),
    ]


# --- Index merge safety: RDKit/FPSim2 upgrades and failed merges -------------

OLD_RDKIT = "2000.01.1"


def _build_old_chunk(paths, sdf, monkeypatch):
    """Build and finalize one chunk as if an older RDKit had fingerprinted it."""
    import rdkit

    with monkeypatch.context() as m:
        m.setattr(rdkit, "__version__", OLD_RDKIT)
        build_chunks(paths, sources([sdf]), workers=1)
        finalize(paths)
    return paths.load_manifest()


def _index_files(paths):
    return sorted(p.name for p in paths.index_dir.iterdir())


def test_chunk_stats_record_versions(built_home):
    from importlib.metadata import version

    import rdkit

    stats = IndexPaths(built_home).completed_chunks()[0]
    assert stats.rdkit_version == rdkit.__version__
    assert stats.fpsim2_version == version("fpsim2")


def test_extending_after_an_rdkit_upgrade_is_refused_and_keeps_the_index(
    tmp_path, sample_sdfs, monkeypatch
):
    paths = IndexPaths(tmp_path)
    before = _build_old_chunk(paths, sample_sdfs[0], monkeypatch)
    with pytest.raises(IndexBuildError, match=f"RDKit {OLD_RDKIT}") as err:
        build_chunks(paths, sources([sample_sdfs[1]]), workers=1)
    assert "still searchable" in str(err.value)
    assert not paths.chunk_stats(CHUNK_NAMES[1]).exists()
    assert paths.load_manifest() == before
    found = service.similarity_search(ASPIRIN, threshold=0.3, home=tmp_path)
    assert found.hits[0].cid == 2244


def _h5_payload(path):
    """Every array and table in an FPSim2 file, as bytes, keyed by node path.
    The config node holds in-memory pointers, so it is left out."""
    with tb.open_file(str(path), mode="r") as f:
        return {
            node._v_pathname: np.asarray(node.read()).tobytes()
            for node in f.walk_nodes("/")
            if hasattr(node, "read") and node._v_name != "config"
        }


def test_unchanged_chunk_from_another_rdkit_is_rebuilt_when_extending(
    tmp_path, sample_sdfs, monkeypatch
):
    """An unchanged local file is rebuilt from its SDF, not refused, after an
    RDKit upgrade: the index can grow without a new data directory."""
    import rdkit

    paths = IndexPaths(tmp_path / "grown")
    _build_old_chunk(paths, sample_sdfs[0], monkeypatch)
    build_chunks(paths, sources(sample_sdfs[:2]), workers=1)
    finalize(paths)

    chunks = {c.name: c for c in paths.completed_chunks()}
    assert {c.rdkit_version for c in chunks.values()} == {rdkit.__version__}
    assert len(chunks) == 2
    fresh = IndexPaths(tmp_path / "fresh")
    build_chunks(fresh, sources(sample_sdfs[:1]), workers=1)
    finalize(fresh)
    expected = fresh.completed_chunks()[0]
    rebuilt = chunks[CHUNK_NAMES[0]]
    assert rebuilt.n_records == expected.n_records == len(CHUNK_A)
    assert rebuilt.n_failed == expected.n_failed
    assert _h5_payload(paths.chunk_h5(CHUNK_NAMES[0])) == _h5_payload(
        fresh.chunk_h5(CHUNK_NAMES[0])
    )
    assert h5_versions(paths.chunk_h5(CHUNK_NAMES[0])) == runtime_versions()
    found = service.similarity_search(ASPIRIN, threshold=0.3, home=tmp_path / "grown")
    assert found.hits[0].cid == 2244


def test_pubchem_chunk_from_another_rdkit_is_still_refused_when_extending(
    tmp_path, sample_sdfs, monkeypatch
):
    """A PubChem chunk has no local file to rebuild from, so it still blocks."""
    paths = IndexPaths(tmp_path)
    before = _build_old_chunk(paths, sample_sdfs[0], monkeypatch)
    pubchem = [SdfSource(name=CHUNK_NAMES[0])]
    with pytest.raises(IndexBuildError, match=f"RDKit {OLD_RDKIT}"):
        build_chunks(paths, pubchem + sources([sample_sdfs[1]]), workers=1)
    assert not paths.chunk_stats(CHUNK_NAMES[1]).exists()
    assert paths.load_manifest() == before


def test_chunk_whose_file_is_gone_from_another_rdkit_is_still_refused(
    tmp_path, sample_sdfs, monkeypatch
):
    """A chunk whose recorded file was deleted cannot be rebuilt, so it blocks."""
    gone_dir = tmp_path / "gone"
    gone_dir.mkdir()
    gone = gone_dir / sample_sdfs[0].name
    gone.write_bytes(sample_sdfs[0].read_bytes())
    paths = IndexPaths(tmp_path / "home")
    before = _build_old_chunk(paths, gone, monkeypatch)
    gone.unlink()
    with pytest.raises(IndexBuildError, match=f"RDKit {OLD_RDKIT}") as err:
        build_chunks(paths, sources([sample_sdfs[1]]), workers=1)
    assert not paths.chunk_stats(CHUNK_NAMES[1]).exists()
    assert paths.load_manifest() == before
    # The message names the chunk, says which file is needed, and keeps the
    # advice that works without the file.
    message = str(err.value)
    assert "index build --sdf FILE" in message
    assert f"with the file of {CHUNK_NAMES[0]} in place of FILE" in message
    assert f"The index recorded {gone.resolve()}, not there now" in message
    assert "PubChem's download" not in message  # it records a file of the user's
    assert "index build --home NEW_DIR" in message


def test_version_refusal_names_the_build_that_brings_the_chunks_up_to_date(
    tmp_path, sample_sdfs, monkeypatch
):
    """Issue #214: the refusal offered only a new data directory or deleting the
    index, although building from the chunks' own files again is enough."""
    paths = IndexPaths(tmp_path)
    _build_old_chunk(paths, sample_sdfs[0], monkeypatch)
    with pytest.raises(IndexBuildError, match="still searchable") as err:
        build_chunks(paths, sources([sample_sdfs[1]]), workers=1)
    message = str(err.value)
    assert f"index build {_sdf_options(sample_sdfs[:1])}" in message
    assert "Then run this build again" in message
    assert "NEW_DIR" not in message

    # The advice works: that build, then the one that was refused.
    build_chunks(paths, sources([sample_sdfs[0]]), workers=1)
    finalize(paths)
    build_chunks(paths, sources([sample_sdfs[1]]), workers=1)
    assert finalize(paths).n_compounds == N_GOOD


def _sdf_options(files):
    """The --sdf options a message names, quoted as the message quotes them."""
    from pubchem_finder import hints

    return " ".join(f"--sdf {hints._quote(str(Path(f).resolve()))}" for f in files)


def test_version_refusal_names_every_stale_chunk_also_those_being_rebuilt(
    tmp_path, sample_sdfs, monkeypatch
):
    """The build that is named has to work by itself. Naming only the chunks the
    refused build left alone gave a build that was refused in turn, for the
    chunks it left alone, and the two refusals pointed at each other."""
    import rdkit

    paths = IndexPaths(tmp_path / "home")
    libs = []
    for name, source in (("liba", sample_sdfs[0]), ("libb", sample_sdfs[1])):
        lib = tmp_path / f"{name}.sdf.gz"
        lib.write_bytes(source.read_bytes())
        libs.append(lib)
    with monkeypatch.context() as m:
        m.setattr(rdkit, "__version__", OLD_RDKIT)
        build_chunks(paths, sources(libs), workers=1)
        finalize(paths)
    new = tmp_path / "libc.sdf.gz"
    new.write_bytes(sample_sdfs[1].read_bytes())

    # Any build that names some of the stale files, but not all, gets the same
    # answer: both of them.
    for asked in ([libs[0]], [libs[1]], [libs[0], new], [new]):
        with pytest.raises(IndexBuildError, match=r"^2 chunk\(s\)") as err:
            build_chunks(paths, sources(asked), workers=1)
        assert f"index build {_sdf_options(libs)}" in str(err.value)
    assert {c.rdkit_version for c in paths.completed_chunks()} == {OLD_RDKIT}

    report = build_chunks(paths, sources(libs), workers=1)
    assert sorted(report.rebuilt) == ["liba", "libb"]
    assert {c.rdkit_version for c in paths.completed_chunks()} == {rdkit.__version__}
    build_chunks(paths, sources([libs[0], new]), workers=1)
    assert sorted(finalize(paths).chunks) == ["liba", "libb", "libc"]


def test_stale_chunk_whose_file_moved_is_rebuilt_from_where_it_is_now(
    tmp_path, sample_sdfs, monkeypatch
):
    """A file that moved was taken for the chunk's own (same size and time), yet
    its stale chunk blocked the build that named it."""
    import rdkit

    first = tmp_path / "data" / "lib.sdf.gz"
    first.parent.mkdir()
    first.write_bytes(sample_sdfs[0].read_bytes())
    paths = IndexPaths(tmp_path / "home")
    _build_old_chunk(paths, first, monkeypatch)
    moved = tmp_path / "data2" / "lib.sdf.gz"
    moved.parent.mkdir()
    first.rename(moved)

    report = build_chunks(paths, sources([moved, sample_sdfs[1]]), workers=1)
    assert report.rebuilt == ["lib"]
    chunks = {c.name: c for c in paths.completed_chunks()}
    assert chunks["lib"].rdkit_version == rdkit.__version__
    assert chunks["lib"].source is not None
    assert chunks["lib"].source.path == str(moved.resolve())
    assert finalize(paths).n_compounds == N_GOOD


def test_version_refusal_asks_for_the_file_of_a_chunk_that_records_none(
    tmp_path, sample_sdfs, monkeypatch
):
    """A chunk built by 2.2.0 or earlier records no file. It was called a PubChem
    chunk, or one whose file is gone, while giving its file rebuilds it."""
    lib = tmp_path / "mylib.sdf.gz"
    lib.write_bytes(sample_sdfs[0].read_bytes())
    paths = IndexPaths(tmp_path / "home")
    _build_old_chunk(paths, lib, monkeypatch)
    stats_file = paths.chunk_stats("mylib")
    stats = json.loads(stats_file.read_text())
    del stats["source"]
    stats_file.write_text(json.dumps(stats))

    with pytest.raises(IndexBuildError, match="still searchable") as err:
        build_chunks(paths, sources([sample_sdfs[1]]), workers=1)
    message = str(err.value)
    assert "index build --sdf FILE" in message
    assert "with the file of mylib in place of FILE" in message
    assert "PubChem's download" not in message
    build_chunks(paths, sources([lib]), workers=1)
    build_chunks(paths, sources([sample_sdfs[1]]), workers=1)
    assert finalize(paths).n_compounds == N_GOOD


def test_version_refusal_for_downloaded_chunks_depends_on_their_sdf_files(
    tmp_path, sample_sdfs, monkeypatch
):
    """A chunk from PubChem's download records no file and is not downloaded
    again. Without its SDF the new directory is the way, said once however many
    such chunks there are (359 for all of PubChem). With the SDF that --keep-sdf
    left, it is fingerprinted again like any other file."""
    import rdkit

    paths = IndexPaths(tmp_path)
    with monkeypatch.context() as m:
        m.setattr(rdkit, "__version__", OLD_RDKIT)
        build_chunks(paths, sources(sample_sdfs[:2]), workers=1)
        finalize(paths)
    for name in CHUNK_NAMES[:2]:
        stats_file = paths.chunk_stats(name)
        stats = json.loads(stats_file.read_text())
        del stats["source"]  # a chunk PubChem's download built records no file
        stats_file.write_text(json.dumps(stats))
    new = tmp_path / "mylib.sdf.gz"
    new.write_bytes(sample_sdfs[1].read_bytes())

    with pytest.raises(IndexBuildError, match="still searchable") as err:
        build_chunks(paths, sources([new]), workers=1)
    message = str(err.value)
    assert (
        f"2 of them came from PubChem's download ({CHUNK_NAMES[0]} to "
        f"{CHUNK_NAMES[1]}) and have no SDF file in {paths.downloads}"
    ) in message
    assert "index build --home NEW_DIR" in message
    assert "--sdf FILE" not in message

    import gzip

    paths.downloads.mkdir()
    kept = [
        paths.downloads / sample_sdfs[0].name,  # as it was downloaded
        paths.downloads / f"{CHUNK_NAMES[1]}.sdf",  # and unpacked
    ]
    kept[0].write_bytes(sample_sdfs[0].read_bytes())
    kept[1].write_bytes(gzip.decompress(sample_sdfs[1].read_bytes()))
    with pytest.raises(IndexBuildError, match="still searchable") as err:
        build_chunks(paths, sources([new]), workers=1)
    message = str(err.value)
    assert f"index build {_sdf_options(kept)}" in message
    assert "PubChem's download" not in message
    build_chunks(paths, sources(kept), workers=1)
    assert {c.rdkit_version for c in paths.completed_chunks()} == {rdkit.__version__}


def _stale_libraries(tmp_path, sample_sdfs, monkeypatch, names):
    """A home whose chunks ``names`` another RDKit fingerprinted, and their files."""
    import rdkit

    data = tmp_path / "data"
    data.mkdir()
    libs = []
    for n, name in enumerate(names):
        lib = data / f"{name}.sdf.gz"
        lib.write_bytes(sample_sdfs[n % 2].read_bytes())
        libs.append(lib)
    paths = IndexPaths(tmp_path / "home")
    with monkeypatch.context() as m:
        m.setattr(rdkit, "__version__", OLD_RDKIT)
        build_chunks(paths, sources(libs), workers=1)
        finalize(paths)
    return paths, libs


def test_version_refusal_names_a_given_file_where_it_is_given(
    tmp_path, sample_sdfs, monkeypatch
):
    """The files of the refused build count: one that moved is named where it is
    now, not asked for, and one given by a relative path is named in full, so
    that the build runs from any directory."""
    paths, (liba, libb) = _stale_libraries(tmp_path, sample_sdfs, monkeypatch, "ab")
    moved = tmp_path / "data2" / liba.name
    moved.parent.mkdir()
    liba.rename(moved)
    monkeypatch.chdir(moved.parent)
    with pytest.raises(IndexBuildError, match="still searchable") as err:
        build_chunks(paths, sources([Path(moved.name)]), workers=1)
    message = str(err.value)
    assert f"index build {_sdf_options([moved, libb])}" in message
    assert "--sdf FILE" not in message
    assert "Then run this build again." in message


def test_version_refusal_names_every_chunk_that_needs_its_file(
    tmp_path, sample_sdfs, monkeypatch
):
    """A build that leaves one out is refused again: no name is cut from the list."""
    paths, libs = _stale_libraries(tmp_path, sample_sdfs, monkeypatch, "abcde")
    for lib in libs[1:]:
        lib.unlink()
    with pytest.raises(IndexBuildError, match=r"^5 chunk\(s\)") as err:
        build_chunks(paths, sources([sample_sdfs[1]]), workers=1)
    message = str(err.value)
    assert f"index build {_sdf_options(libs[:1])} --sdf FILE --sdf FILE" in message
    assert message.count("--sdf FILE") == 4
    assert "with the files of b, c, d, e in place of each FILE" in message
    assert "(b.sdf or b.sdf.gz for b)" in message
    recorded = " and ".join(str(lib.resolve()) for lib in libs[1:])
    assert f"The index recorded {recorded}, not there now" in message
    assert "under another name" not in message  # said of files that are there


def test_version_refusal_names_each_version_the_chunks_were_built_with(
    tmp_path, sample_sdfs, monkeypatch
):
    import rdkit

    paths, _ = _stale_libraries(tmp_path, sample_sdfs, monkeypatch, "a")
    older = tmp_path / "data" / "b.sdf.gz"
    older.write_bytes(sample_sdfs[1].read_bytes())
    with monkeypatch.context() as m:
        m.setattr(rdkit, "__version__", "1999.09.1")
        m.setattr(builder_mod, "_check_versions", lambda *args, **kwargs: None)
        build_chunks(paths, sources([older]), workers=1)
    with pytest.raises(IndexBuildError) as err:
        build_chunks(paths, sources([sample_sdfs[1]]), workers=1)
    assert (
        f"2 chunk(s) in {paths.chunks} were fingerprinted with RDKit {OLD_RDKIT} / "
        f"FPSim2 {runtime_versions()[1]} and RDKit 1999.09.1 / FPSim2 "
    ) in str(err.value)


def test_version_refusal_names_the_file_of_a_chunk_named_by_an_older_release(
    tmp_path, sample_sdfs, monkeypatch
):
    """Releases before 2.2.2 kept an upper-case suffix in the chunk's name
    (`Library.SDF` for Library.SDF.gz). Such a chunk's file is still its file."""
    import rdkit

    sdf = tmp_path / "Library.SDF.GZ"
    sdf.write_bytes(sample_sdfs[0].read_bytes())
    paths = IndexPaths(tmp_path / "home")
    with monkeypatch.context() as m:
        m.setattr(rdkit, "__version__", OLD_RDKIT)
        build_chunks(
            paths, [SdfSource(name="Library.SDF.GZ", local_path=sdf)], workers=1
        )
        finalize(paths)
    assert [c.name for c in paths.completed_chunks()] == ["Library.SDF.GZ"]
    with pytest.raises(IndexBuildError, match="still searchable") as err:
        build_chunks(paths, sources([sample_sdfs[1]]), workers=1)
    message = str(err.value)
    assert f"index build {_sdf_options([sdf])}" in message
    assert "--sdf FILE" not in message
    build_chunks(paths, sources([sdf]), workers=1)
    assert {c.rdkit_version for c in paths.completed_chunks()} == {rdkit.__version__}


def test_version_refusal_names_a_file_only_under_its_chunks_name(
    tmp_path, sample_sdfs, monkeypatch
):
    """A file is taken for its chunk by its name. One given through a link of
    another name (a `latest` link, a checkout that names files by their hash) is
    recorded under the name it resolves to: named in the build, it would become
    a chunk of its own and the refusal would come again, word for word."""
    data = tmp_path / "data"
    data.mkdir()
    target = data / "vendor_2026_11.sdf.gz"
    target.write_bytes(sample_sdfs[0].read_bytes())
    link = data / "current.sdf.gz"
    try:
        link.symlink_to(target)
    except OSError:
        pytest.skip("no permission to create links here")
    paths = IndexPaths(tmp_path / "home")
    _build_old_chunk(paths, link, monkeypatch)
    assert [c.name for c in paths.completed_chunks()] == ["current"]

    # Not given: its recorded file is not named in the build, but asked for.
    with pytest.raises(IndexBuildError, match="still searchable") as err:
        build_chunks(paths, sources([sample_sdfs[1]]), workers=1)
    message = str(err.value)
    assert "index build --sdf FILE" in message
    assert "vendor_2026_11.sdf.gz --home" not in message
    assert "with the file of current in place of FILE" in message
    assert f"The index recorded {target.resolve()}, under another name" in message
    assert "not there now" not in message  # it is there

    # Given: it is rebuilt, and with another stale chunk beside it the build
    # that is named has the link in it, in full.
    other = data / "libb.sdf.gz"
    other.write_bytes(sample_sdfs[1].read_bytes())
    _build_old_chunk(paths, other, monkeypatch)
    monkeypatch.chdir(data)
    with pytest.raises(IndexBuildError, match="still searchable") as err:
        build_chunks(paths, sources([Path(link.name)]), workers=1)
    quoted = " ".join(f"--sdf {hints._quote(str(p))}" for p in (link, other))
    assert f"index build {quoted}" in str(err.value)

    report = build_chunks(paths, sources([link, other]), workers=1)
    assert sorted(report.rebuilt) == ["current", "libb"]


@pytest.mark.skipif(
    os.name == "nt", reason="Windows takes `..` off a path before it follows links"
)
def test_version_refusal_names_a_file_given_through_a_linked_folder(
    tmp_path, sample_sdfs, monkeypatch
):
    """`..` after a folder that is a link is the parent of the folder it links
    to. That is where the system finds the file, not what is left when the two
    names are struck from the path, which the refusal named."""
    paths, (liba, libb) = _stale_libraries(tmp_path, sample_sdfs, monkeypatch, "ab")
    elsewhere = tmp_path / "elsewhere"
    (elsewhere / "sub").mkdir(parents=True)
    (tmp_path / "proj").mkdir()
    (tmp_path / "proj" / "sub").symlink_to(elsewhere / "sub")
    moved = elsewhere / libb.name
    libb.rename(moved)
    monkeypatch.chdir(tmp_path / "proj")
    through = Path("sub") / ".." / libb.name
    new = tmp_path / "new.sdf.gz"
    new.write_bytes(sample_sdfs[1].read_bytes())
    with pytest.raises(IndexBuildError, match="still searchable") as err:
        build_chunks(paths, sources([through, new]), workers=1)
    named = (liba, moved)
    quoted = " ".join(f"--sdf {hints._quote(str(p))}" for p in named)
    assert f"index build {quoted}" in str(err.value)
    assert "--sdf FILE" not in str(err.value)
    report = build_chunks(paths, sources(named), workers=1)
    assert sorted(report.rebuilt) == ["a", "b"]


@pytest.mark.skipif(os.name == "nt", reason="POSIX directory modes")
def test_version_refusal_survives_a_recorded_file_that_cannot_be_reached(
    tmp_path, sample_sdfs, monkeypatch
):
    """Before Python 3.14 pathlib raises for a file in a directory that may not
    be entered, and the refusal came out as a bare "Permission denied"."""
    paths, (liba,) = _stale_libraries(tmp_path, sample_sdfs, monkeypatch, "a")
    probe = tmp_path / "probe"
    probe.mkdir()
    probe.chmod(0)
    held_to_modes = not os.access(probe, os.R_OK | os.X_OK)
    probe.chmod(0o755)
    if not held_to_modes:
        pytest.skip("this user is not held to file modes")
    liba.parent.chmod(0)
    try:
        with pytest.raises(IndexBuildError, match="still searchable") as err:
            build_chunks(paths, sources([sample_sdfs[1]]), workers=1)
    finally:
        liba.parent.chmod(0o755)
    assert "with the file of a in place of FILE" in str(err.value)


def test_chunks_from_2_1_0_are_checked_through_their_h5(
    tmp_path, sample_sdfs, monkeypatch
):
    """Chunks built before 2.1.1 record no versions in their JSON."""
    paths = IndexPaths(tmp_path)
    _build_old_chunk(paths, sample_sdfs[0], monkeypatch)
    stats_file = paths.chunk_stats(CHUNK_NAMES[0])
    stats = json.loads(stats_file.read_text())
    del stats["rdkit_version"], stats["fpsim2_version"]
    stats_file.write_text(json.dumps(stats))
    with pytest.raises(IndexBuildError, match=f"RDKit {OLD_RDKIT}"):
        build_chunks(paths, sources([sample_sdfs[1]]), workers=1)


def test_finalize_refuses_mixed_versions_before_touching_the_index(
    tmp_path, sample_sdfs, monkeypatch
):
    from pubchem_finder.index.builder import build_chunk

    paths = IndexPaths(tmp_path)
    before = _build_old_chunk(paths, sample_sdfs[0], monkeypatch)
    fingerprints = paths.fingerprints.read_bytes()
    build_chunk(str(sample_sdfs[1]), str(tmp_path))  # skips build_chunks' check
    with pytest.raises(IndexBuildError, match="different versions"):
        finalize(paths)
    assert paths.load_manifest() == before
    assert paths.fingerprints.read_bytes() == fingerprints
    assert _index_files(paths) == ["manifest.json", "pubchem.h5"]


def test_failed_merge_keeps_the_previous_index(tmp_path, sample_sdfs, monkeypatch):
    from pubchem_finder.index import builder

    paths = IndexPaths(tmp_path)
    build_chunks(paths, sources([sample_sdfs[0]]), workers=1)
    before = finalize(paths)
    build_chunks(paths, sources([sample_sdfs[1]]), workers=1)

    def full_disk(src, dest):
        raise OSError(28, "No space left on device")

    monkeypatch.setattr(builder, "sort_fpsim2_file", full_disk)
    with pytest.raises(IndexBuildError, match=r"No space left.*unchanged"):
        finalize(paths)
    assert paths.load_manifest() == before
    assert _index_files(paths) == ["manifest.json", "pubchem.h5"]


def test_failed_swap_restores_the_manifest(tmp_path, sample_sdfs, monkeypatch):
    """On Windows a running server can hold pubchem.h5 open, so replacing fails."""
    import os

    paths = IndexPaths(tmp_path)
    build_chunks(paths, sources([sample_sdfs[0]]), workers=1)
    before = finalize(paths)
    build_chunks(paths, sources([sample_sdfs[1]]), workers=1)
    real_replace = os.replace

    def locked(src, dst):
        if str(dst).endswith("pubchem.h5"):
            raise PermissionError(13, "The process cannot access the file")
        real_replace(src, dst)

    monkeypatch.setattr(os, "replace", locked)
    with pytest.raises(IndexBuildError, match="stop it and run the build again"):
        finalize(paths)
    monkeypatch.undo()
    assert paths.load_manifest() == before
    assert _index_files(paths) == ["manifest.json", "pubchem.h5"]


def test_unreadable_sdf_is_an_index_build_error(tmp_path, sample_sdfs):
    data = sample_sdfs[1].read_bytes()
    broken = tmp_path / f"{CHUNK_NAMES[1]}.sdf.gz"
    broken.write_bytes(data[: len(data) // 2])
    with pytest.raises(IndexBuildError) as err:
        build_chunks(IndexPaths(tmp_path / "home"), sources([broken]), workers=1)
    assert f"{broken} is damaged or incomplete" in str(err.value)


# --- Overlapping chunks: each CID once ----------------------------------------


@pytest.fixture(scope="module")
def overlap_home(tmp_path_factory, sample_sdfs):
    """Chunk A plus an SDF that repeats three of its compounds."""
    from tests.conftest import CHUNK_A, write_sdf

    d = tmp_path_factory.mktemp("overlap")
    mine = write_sdf(d / "my_compounds.sdf", CHUNK_A[:3], gz=False)
    paths = IndexPaths(d / "home")
    build_chunks(paths, sources([sample_sdfs[0], mine]), workers=1)
    finalize(paths, substructure=True)
    return paths.home


def test_overlapping_chunks_give_each_cid_once(overlap_home):
    found = service.similarity_search(
        ASPIRIN, threshold=0.01, top_k=1000, home=overlap_home
    )
    cids = [h.cid for h in found.hits]
    assert cids[0] == 2244
    assert len(cids) == len(set(cids)) == 5
    assert found.truncated is False
    top3 = service.similarity_search(
        ASPIRIN, threshold=0.01, top_k=3, home=overlap_home
    )
    assert len({h.cid for h in top3.hits}) == 3
    assert top3.truncated is True
    # Eight rows reach the threshold, for five compounds: a list of five is whole.
    top5 = service.similarity_search(
        ASPIRIN, threshold=0.01, top_k=5, home=overlap_home
    )
    assert len({h.cid for h in top5.hits}) == 5
    assert top5.truncated is False
    assert top5.note is None


@pytest.mark.parametrize(
    "raised",
    ["RuntimeError: Query interrupted", "duckdb.InterruptException", "another error"],
)
def test_ctrl_c_during_the_metadata_query_is_ctrl_c(built_home, raised):
    """DuckDB turns a Ctrl+C that arrives while it runs a query into an error of
    its own, or a bare RuntimeError. Taken for a failed query, it made a batch
    search throw away the rows it keeps on Ctrl+C, about one time in three."""
    import duckdb

    index = LocalIndex(built_home)
    chunks = list(index.merged.manifest.chunks)
    errors = {
        "RuntimeError: Query interrupted": RuntimeError("Query interrupted"),
        "duckdb.InterruptException": duckdb.InterruptException("Interrupted!"),
        "another error": RuntimeError("out of memory"),
    }

    class Interrupted:
        def cursor(self):
            return self

        def execute(self, sql):
            raise errors[raised]

    index._db = Interrupted()
    index._db_files = index.paths.metadata_files(chunks)
    expected = RuntimeError if raised == "another error" else KeyboardInterrupt
    with pytest.raises(expected):
        index._metadata([2244], chunks)
    if raised != "another error":
        # And through the search, which tries a failed query again: not this one.
        with pytest.raises(KeyboardInterrupt):
            index.metadata([2244])


def test_a_failed_query_that_mentions_interruption_is_not_ctrl_c(built_home):
    """Told by what the error is, not by a word in it: DuckDB's errors quote
    paths, and a data directory may be called anything. Such an error is an
    error of the query, which is tried again; as Ctrl+C it ended an MCP server."""
    import duckdb

    index = LocalIndex(built_home)
    chunks = list(index.merged.manifest.chunks)
    failed = duckdb.IOException(
        'Cannot open file "/data/interrupted-screens/home/chunks/lib1.parquet"'
    )
    for words in (
        'Cannot open file "/data/Query interrupted/x.parquet"',
        "Query interrupted by a damaged footer",
    ):
        index = LocalIndex(built_home)

        class Wordy:
            def cursor(self):
                return self

            def execute(self, sql, words=words):
                raise RuntimeError(words)

        index._db = Wordy()
        index._db_files = index.paths.metadata_files(list(index.merged.manifest.chunks))
        caught: BaseException | None = None
        try:
            index._metadata([2244], list(index.merged.manifest.chunks))
        except KeyboardInterrupt:
            pytest.fail(f"taken for Ctrl+C: {words}")
        except RuntimeError as exc:
            caught = exc
        assert str(caught) == words
    calls = []

    class Failing:
        def cursor(self):
            return self

        def execute(self, sql):
            calls.append(sql)
            raise failed

    index._db = Failing()
    index._db_files = index.paths.metadata_files(chunks)
    raised: BaseException | None = None
    try:
        index._metadata([2244], chunks)
    except KeyboardInterrupt:  # caught here: let go, it would end the test run
        pytest.fail("a failed query was taken for Ctrl+C")
    except duckdb.IOException as exc:
        raised = exc
    assert raised is failed
    # The search tries again, with a connection of its own, and answers.
    index._db = Failing()
    try:
        found = index.metadata([2244])
    except KeyboardInterrupt:
        pytest.fail("a failed query was taken for Ctrl+C")
    assert 2244 in found
    assert len(calls) == 2


def test_overlapping_chunks_give_each_substructure_hit_once(overlap_home):
    found = service.substructure_search("c1ccccc1C(=O)O", home=overlap_home)
    assert [h.cid for h in found.hits] == [243, 338, 2244, 4133]


def test_finalize_refuses_duplicate_id_when_invalid_smiles_tag_is_replaced(tmp_path):
    """Issue #169: an invalid tag must not hide a different structure sharing its ID."""
    from rdkit import Chem

    first_path = tmp_path / "tagged.sdf"
    second_path = tmp_path / "untagged.sdf"
    for path, smiles, invalid_tag in (
        (first_path, "CCO", True),
        (second_path, "CCC", False),
    ):
        molecule = Chem.MolFromSmiles(smiles)
        molecule.SetProp("PUBCHEM_COMPOUND_CID", "4242")
        if invalid_tag:
            molecule.SetProp("PUBCHEM_SMILES", "not a smiles")
        with Chem.SDWriter(str(path)) as writer:
            writer.write(molecule)

    paths = IndexPaths(tmp_path / "home")
    build_chunks(paths, sources([first_path, second_path]), workers=1)

    with pytest.raises(IndexBuildError, match="ID 4242 stands for different compounds"):
        finalize(paths)


def test_search_asks_for_more_when_duplicates_crowd_out_hits(built_home, monkeypatch):
    import numpy as np
    from rdkit import Chem

    rows = np.array(
        [
            (2244, 1.0),
            (2244, 1.0),
            (338, 0.8),
            (338, 0.8),
            (4133, 0.7),
            (4133, 0.7),
            (243, 0.6),
        ],
        dtype=[("mol_id", "<u4"), ("coeff", "<f4")],
    )
    asked: list[int] = []

    class StubEngine:
        def top_k(self, query, k, **kwargs):
            asked.append(k)
            return rows[:k]

    index = LocalIndex(built_home, mode="memory")
    monkeypatch.setattr(index, "_ensure_loaded", lambda: (StubEngine(), index.merged))
    index._in_memory = True
    assert index.search(Chem.MolFromSmiles(ASPIRIN), 0.5, 3) == [
        (2244, 1.0),
        (338, 0.8),
        (4133, 0.7),
    ]
    # The first ask is one past the cut; the duplicate 2244/338 copies leave only
    # two distinct CIDs in hand, so the search asks again and exhausts the hits.
    assert asked == [4, 8]


def test_search_widens_once_past_a_tie_at_the_cut(built_home, monkeypatch):
    """A tie at the cut is complete once the lowest score in hand is below it.

    The first fetch ends on a tie, so the score after the cut is no proof that
    nothing outside the fetch ties with it; the second fetch shows a lower score,
    which is. The search must stop there: two calls, not a doubling all the way
    to the end of the hits (#159)."""
    import numpy as np
    from rdkit import Chem

    rows = np.array(
        [
            (450661, 1.0),
            (2244, 1.0),
            (243, 0.8),
            (338, 0.7),
        ],
        dtype=[("mol_id", "<u4"), ("coeff", "<f4")],
    )
    asked: list[int] = []

    class StubEngine:
        def top_k(self, query, k, **kwargs):
            asked.append(k)
            return rows[:k]

    index = LocalIndex(built_home, mode="memory")
    monkeypatch.setattr(index, "_ensure_loaded", lambda: (StubEngine(), index.merged))
    index._in_memory = True
    assert index.search(Chem.MolFromSmiles(ASPIRIN), 0.5, 1) == [(2244, 1.0)]
    assert asked == [2, 4]


# --- Index format ------------------------------------------------------------


def _copy_home(built_home, tmp_path):
    import shutil

    home = tmp_path / "home"
    shutil.copytree(built_home, home)
    return home


def test_manifest_records_the_format_version(built_home):
    from pubchem_finder.index.layout import FORMAT_VERSION

    manifest = json.loads(IndexPaths(built_home).manifest.read_text())
    assert manifest["format_version"] == FORMAT_VERSION


def test_manifest_from_2_1_reads_as_format_1(tmp_path, built_home):
    home = _copy_home(built_home, tmp_path)
    paths = IndexPaths(home)
    manifest = json.loads(paths.manifest.read_text())
    del manifest["format_version"]
    paths.manifest.write_text(json.dumps(manifest))
    loaded = paths.load_manifest()
    assert loaded is not None
    assert loaded.format_version == 1
    found = service.similarity_search(ASPIRIN, threshold=0.3, home=home)
    assert found.hits[0].cid == 2244


def test_index_in_a_newer_format_is_refused(tmp_path, built_home):
    from pubchem_finder.index.layout import FORMAT_VERSION

    home = _copy_home(built_home, tmp_path)
    paths = IndexPaths(home)
    manifest = json.loads(paths.manifest.read_text())
    manifest["format_version"] = FORMAT_VERSION + 1
    manifest["pubchem_finder_version"] = "9.0.0"
    manifest["fp_params"] = "a shape this version can't validate"
    paths.manifest.write_text(json.dumps(manifest))
    with pytest.raises(IndexIncompatible, match=r"pubchem-finder 9\.0\.0") as err:
        service.similarity_search(ASPIRIN, home=home)
    assert "Upgrade pubchem-finder" in str(err.value)
    with pytest.raises(IndexIncompatible):
        service.index_status(home)


def test_metadata_reads_chunks_with_different_columns(tmp_path, built_home):
    """A later version may add Parquet columns; older chunks lack them."""
    import pyarrow.parquet as pq

    home = _copy_home(built_home, tmp_path)
    path = IndexPaths(home).chunk_parquet(CHUNK_NAMES[1])
    table = pq.read_table(path)
    pq.write_table(table.append_column("title", [["x"] * table.num_rows]), path)
    rows = LocalIndex(home).metadata([2244, 2519])
    assert sorted(rows) == [2244, 2519]
    assert rows[2519]["iupac_name"] == "1,3,7-trimethylpurine-2,6-dione"


# --- Chunks from any reader -------------------------------------------------------


def test_write_chunk_indexes_records_from_any_reader(tmp_path):
    """The demo set comes from a TSV, not an SDF: write_chunk takes any reader."""
    from rdkit import Chem

    from pubchem_finder.data.sdf import SdfRecord
    from pubchem_finder.index.builder import write_chunk

    class Records:
        n_ok = 0
        n_failed = 1

        def __iter__(self):
            for cid, smiles in [(2244, ASPIRIN), (702, "CCO")]:
                self.n_ok += 1
                yield SdfRecord(
                    cid, Chem.MolFromSmiles(smiles), smiles, None, None, None, None
                )

    paths = IndexPaths(tmp_path)
    stats = write_chunk(paths, "drugs", Records())
    assert (stats.name, stats.n_records, stats.n_failed) == ("drugs", 2, 1)
    finalize(paths)
    found = service.similarity_search(ASPIRIN, home=tmp_path)
    assert [h.cid for h in found.hits] == [2244]


# --- FPSim2 threads ------------------------------------------------------------


def test_in_memory_searches_run_in_one_thread(substructure_home, monkeypatch):
    """FPSim2 0.7.4's threaded search updates the fingerprint arrays' reference
    counts unsafely and can free them mid-search; with n_workers=1 it searches in
    the calling thread."""
    import sys

    from rdkit import Chem

    index = LocalIndex(substructure_home, mode="memory")
    query = Chem.MolFromSmiles(ASPIRIN)
    fragment = Chem.MolFromSmiles("c1ccccc1")
    index.search(query, 0.1, 5)
    index.screen(fragment)
    engines = (index._engine, index._pattern_engine)
    workers: list[int] = []
    for engine, method in zip(engines, ("top_k", "substructure"), strict=True):
        real = getattr(engine, method)

        def spy(*args, _real=real, **kwargs):
            workers.append(kwargs["n_workers"])
            return _real(*args, **kwargs)

        monkeypatch.setattr(engine, method, spy)
    before = [sys.getrefcount(e.fps) for e in engines]
    for _ in range(200):
        index.search(query, 0.1, 5)
        index.screen(fragment)
    assert set(workers) == {1}
    assert [sys.getrefcount(e.fps) for e in engines] == before


def test_in_memory_searches_of_one_array_never_overlap(substructure_home, monkeypatch):
    """FPSim2 drops its reference to the fingerprint array without the GIL, so two
    threads searching one in-memory array race on its reference count even with
    n_workers=1. Parallel MCP tool calls do exactly that; LocalIndex serializes them
    per engine."""
    import threading
    import time
    from concurrent.futures import ThreadPoolExecutor

    from rdkit import Chem

    index = LocalIndex(substructure_home, mode="memory")
    query = Chem.MolFromSmiles(ASPIRIN)
    fragment = Chem.MolFromSmiles("c1ccccc1")
    index.search(query, 0.1, 5)
    index.screen(fragment)
    guard = threading.Lock()
    active = {"top_k": 0, "substructure": 0}
    peak = dict(active)
    engines = (index._engine, index._pattern_engine)
    for engine, method in zip(engines, active, strict=True):
        real = getattr(engine, method)

        def spy(*args, _real=real, _method=method, **kwargs):
            with guard:
                active[_method] += 1
                peak[_method] = max(peak[_method], active[_method])
            try:
                time.sleep(0.01)
                return _real(*args, **kwargs)
            finally:
                with guard:
                    active[_method] -= 1

        monkeypatch.setattr(engine, method, spy)
    with ThreadPoolExecutor(max_workers=8) as pool:
        jobs = [pool.submit(index.search, query, 0.1, 5) for _ in range(8)]
        jobs += [pool.submit(index.screen, fragment) for _ in range(8)]
        for job in jobs:
            job.result()
    assert peak == {"top_k": 1, "substructure": 1}


# --- Chunks rebuilt since the merge ---------------------------------------------


def _chunk(name, build_id, source=True, modified_ns=1):
    from pubchem_finder.index.layout import ChunkStats, SourceFile

    return ChunkStats(
        name=name,
        n_records=1,
        n_failed=0,
        seconds=0.0,
        fp_type="Morgan",
        fp_params={"radius": 2, "fpSize": 1024},
        build_id=build_id,
        source=SourceFile(path=f"/data/{name}.sdf", size=1, mtime_ns=modified_ns)
        if source
        else None,
    )


def test_unmerged_compares_build_ids(built_home):
    from pubchem_finder.index.layout import unmerged

    manifest = IndexPaths(built_home).load_manifest()
    held = manifest.model_copy(
        update={
            "chunks": ["a", "b", "c"],
            "chunk_builds": {"a": "1", "b": "2", "c": "3"},
        }
    )
    on_disk = [_chunk("a", "1"), _chunk("b", "9"), _chunk("new", "5")]
    # b is another build, c has no stats; a chunk that was never merged is not
    # the index's business.
    assert unmerged(held, on_disk) == ["b", "c"]
    assert unmerged(held, [_chunk("a", "1"), _chunk("b", "2"), _chunk("c", "3")]) == []
    # 2.2.1 rebuilds b: its stats lose their id.
    assert unmerged(held, [_chunk("a", "1"), _chunk("b", None), _chunk("c", "3")]) == [
        "b"
    ]


def test_unmerged_under_a_manifest_without_build_ids(built_home):
    """Merged by 2.2.1 or earlier, which wrote the manifest last: a chunk of the
    user's own file with newer stats was built since. PubChem's chunks are never
    rebuilt, and a build id says nothing (2.2.1 can merge a chunk that has one)."""
    from pubchem_finder.index.layout import Merged, unmerged

    manifest = IndexPaths(built_home).load_manifest()
    legacy = manifest.model_copy(
        update={"chunks": ["lib", "old", "pubchem"], "chunk_builds": None}
    )
    on_disk = [
        _chunk("lib", "7"),
        _chunk("old", None),
        _chunk("pubchem", "8", source=False),
    ]
    at_100 = Merged(legacy, 100)  # manifest.json written at 100
    assert unmerged(legacy, on_disk) == []
    assert unmerged(at_100, on_disk, {"lib": 50, "old": 100, "pubchem": 50}) == []
    assert unmerged(at_100, on_disk, {"lib": 101, "old": 50, "pubchem": 50}) == ["lib"]
    # 2.2.1 rebuilt "old" and failed before the merge: no build id, newer stats.
    late = {"lib": 101, "old": 102, "pubchem": 103}
    assert unmerged(at_100, on_disk, late) == ["lib", "old"]

    # The same home after a copy that kept no file times, or after 2.2.1 rewrote
    # the manifest to add substructure search: the chunk still records a file
    # modified after the merge, which no chunk in that merge can have been built
    # from.
    after = int(legacy.built_at.timestamp() * 1e9) + 1
    rebuilt = [_chunk("lib", "7"), _chunk("old", None, modified_ns=after), on_disk[2]]
    assert unmerged(legacy, rebuilt) == ["old"]
    copied = Merged(legacy, after + 10)  # every file as new as the copy
    everything = dict.fromkeys(("lib", "old", "pubchem"), after + 5)
    assert unmerged(copied, rebuilt, everything) == ["old"]
    # Unless the file is dated later than the chunk that was built from it: a
    # clock that is ahead, or a date in the future, in a home that is fine.
    assert unmerged(copied, rebuilt, {**everything, "old": after - 1}) == []
    ids = {"lib": "7", "old": None, "pubchem": "8"}
    held = legacy.model_copy(update={"chunk_builds": ids})
    # Only where the manifest has no ids: with them, the ids decide.
    assert unmerged(held, rebuilt) == []


def _written(path, when_ns):
    os.utime(path, ns=(when_ns, when_ns))


def test_chunks_merged_by_2_2_1_after_this_version_built_them_are_not_stale(
    tmp_path, sample_sdfs
):
    """Every chunk built by this version was reported as indexed again, and lost
    its structures, once 2.2.1 had written the manifest: for instance after
    `index build --with-substructure`, the command 2.2.1's own hint prints."""
    paths = IndexPaths(tmp_path / "home")
    build_chunks(paths, sources(sample_sdfs[:1]), workers=1)
    finalize(paths)
    assert "2.2.1 or earlier" not in (service.index_status(paths.home).hint or "")
    data = json.loads(paths.manifest.read_text())
    del data["chunk_builds"]
    paths.manifest.write_text(json.dumps(data))
    stats = paths.chunk_stats(CHUNK_NAMES[0])
    merged_at = stats.stat().st_mtime_ns + 10**9
    _written(paths.manifest, merged_at)

    status = service.index_status(paths.home)
    assert "merged by pubchem-finder 2.2.1 or earlier" in status.hint
    assert status.unmerged_chunks == []
    found = service.similarity_search(ASPIRIN, threshold=0.9, home=paths.home)
    assert [(h.cid, h.smiles) for h in found.hits] == [(2244, ASPIRIN)]
    assert found.warning is None

    # Indexed again after that merge, by any version, and not merged since.
    _written(stats, merged_at + 10**9)
    assert service.index_status(paths.home).unmerged_chunks == [CHUNK_NAMES[0]]
    found = service.similarity_search(ASPIRIN, threshold=0.9, home=paths.home)
    assert [(h.cid, h.smiles) for h in found.hits] == [(2244, None)]
    assert "was indexed again after the last merge" in found.warning
    _written(stats, merged_at)  # written in the same clock tick: not after
    assert service.index_status(paths.home).unmerged_chunks == []


def test_a_stats_file_is_not_trusted_from_memory_while_it_is_new(tmp_path, sample_sdfs):
    """A running server remembers the stats it has read, by inode number, time and
    size. A rebuild within one tick of a coarse clock, at a reused inode number,
    leaves a file that looks the same: for a few seconds it is read again."""
    paths = IndexPaths(tmp_path / "home")
    build_chunks(paths, sources(sample_sdfs[:1]), workers=1)
    finalize(paths)
    index = LocalIndex(paths.home)
    assert index.unmerged() == {}
    stats = paths.chunk_stats(CHUNK_NAMES[0])
    before = stats.stat()
    text = stats.read_text()
    build_id = json.loads(text)["build_id"]
    with stats.open("r+") as handle:  # in place: the same inode and size
        handle.write(text.replace(build_id, "f" * len(build_id)))
    os.utime(stats, ns=(before.st_atime_ns, before.st_mtime_ns))
    after = stats.stat()
    assert (after.st_ino, after.st_mtime_ns, after.st_size) == (
        before.st_ino,
        before.st_mtime_ns,
        before.st_size,
    )
    assert list(index.unmerged()) == [CHUNK_NAMES[0]]


def test_status_has_no_2_2_1_note_for_an_index_of_pubchem_chunks(tmp_path, sample_sdfs):
    """2.2.1 could only miss a rebuilt --sdf file; PubChem's chunks are never
    rebuilt, so their index is not merged again and nothing is said."""
    paths = IndexPaths(tmp_path / "home")
    name = chunk_name(sample_sdfs[0])
    build_chunks(paths, [SdfSource(name=name, fetch=lambda: sample_sdfs[0])])
    finalize(paths)
    data = json.loads(paths.manifest.read_text())
    del data["chunk_builds"]
    paths.manifest.write_text(json.dumps(data))
    _written(paths.chunk_stats(name), paths.manifest.stat().st_mtime_ns + 10**9)
    status = service.index_status(paths.home)
    assert "2.2.1 or earlier" not in (status.hint or "")
    assert status.unmerged_chunks == []


def test_metadata_reads_again_when_a_file_changes_under_it(built_home, monkeypatch):
    """A build replacing a Parquet file while DuckDB read it surfaced as a raw
    "TProtocolException: Invalid data" from a search."""
    import duckdb

    index = LocalIndex(built_home)
    real = index._metadata
    failures = [duckdb.IOException("Invalid data"), duckdb.IOException("Invalid data")]

    def flaky(cids, chunks):
        if failures:
            raise failures.pop()
        return real(cids, chunks)

    monkeypatch.setattr(index, "_metadata", flaky)
    assert index.metadata([2244])[2244]["smiles"] == ASPIRIN
    failures.extend(duckdb.IOException("Invalid data") for _ in range(3))
    with pytest.raises(duckdb.IOException):
        index.metadata([2244])
