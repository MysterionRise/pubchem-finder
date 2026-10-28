"""Tests for building and searching the local index."""

import json

import pytest

from pubchem_finder import service
from pubchem_finder.core.exceptions import (
    IndexBuildError,
    IndexNotBuilt,
    MoleculeValidationError,
    PubChemFinderError,
)
from pubchem_finder.index import engine as engine_mod
from pubchem_finder.index.builder import SdfSource, build_chunks, finalize
from pubchem_finder.index.engine import LocalIndex, choose_in_memory
from pubchem_finder.index.layout import IndexPaths, chunk_name
from tests.conftest import CHUNK_NAMES, N_GOOD

ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"


def sources(paths):
    return [SdfSource(name=chunk_name(p), local_path=p) for p in paths]


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
    monkeypatch.setenv("PUBCHEM_FINDER_INDEX_MODE", "bogus")
    assert LocalIndex(built_home).mode == "auto"


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
    monkeypatch.setattr(index, "_ensure_loaded", StubEngine)
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
    with pytest.raises(IndexBuildError, match=f"Building {CHUNK_NAMES[1]} failed"):
        build_chunks(IndexPaths(tmp_path / "home"), sources([broken]), workers=1)


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
    top3 = service.similarity_search(
        ASPIRIN, threshold=0.01, top_k=3, home=overlap_home
    )
    assert len({h.cid for h in top3.hits}) == 3


def test_overlapping_chunks_give_each_substructure_hit_once(overlap_home):
    found = service.substructure_search("c1ccccc1C(=O)O", home=overlap_home)
    assert [h.cid for h in found.hits] == [243, 338, 2244, 4133]


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
    monkeypatch.setattr(index, "_ensure_loaded", StubEngine)
    index._in_memory = True
    assert index.search(Chem.MolFromSmiles(ASPIRIN), 0.5, 3) == [
        (2244, 1.0),
        (338, 0.8),
        (4133, 0.7),
    ]
    assert asked == [3, 6]
