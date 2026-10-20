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
from pubchem_finder.index.engine import SimilarityIndex, choose_in_memory
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
    mem = SimilarityIndex(built_home, mode="memory")
    disk = SimilarityIndex(built_home, mode="disk")
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
    assert SimilarityIndex(built_home).metadata([]) == {}


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
    assert choose_in_memory(manifest, "memory") is True
    assert choose_in_memory(manifest, "disk") is False
    monkeypatch.setattr(engine_mod, "available_memory", lambda: 10)
    assert choose_in_memory(manifest, "auto") is False
    monkeypatch.setattr(engine_mod, "available_memory", lambda: None)
    assert choose_in_memory(manifest, "auto") is True


def test_mode_from_environment(monkeypatch, built_home):
    monkeypatch.setenv("PUBCHEM_FINDER_INDEX_MODE", "disk")
    assert SimilarityIndex(built_home).mode == "disk"
    monkeypatch.setenv("PUBCHEM_FINDER_INDEX_MODE", "bogus")
    assert SimilarityIndex(built_home).mode == "auto"


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
