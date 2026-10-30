"""Local substructure search: pattern index build, screen + verify, removal."""

import json
import shutil

import pytest
from rdkit import Chem

from pubchem_finder import service
from pubchem_finder.core.exceptions import (
    IndexBuildError,
    IndexNotBuilt,
    MoleculeValidationError,
    PubChemFinderError,
)
from pubchem_finder.index.builder import (
    SdfSource,
    build_chunks,
    build_pattern_chunks,
    finalize,
)
from pubchem_finder.index.engine import LocalIndex
from pubchem_finder.index.layout import IndexPaths, Manifest, chunk_name
from tests.conftest import CHUNK_A, CHUNK_B, CHUNK_NAMES, N_GOOD

ALL = sorted((cid, smiles) for cid, smiles, *_ in CHUNK_A + CHUNK_B)
INDEXED_CHUNKS = list(CHUNK_NAMES[:2])  # chunk C has no parseable records

QUERIES = [
    ("c1ccccc1", "smiles"),
    ("c1ccccc1C(=O)O", "smiles"),
    ("C(=O)O", "smiles"),
    ("O", "smiles"),
    ("C1CCCCC1", "smiles"),  # no matches
    ("[CX3](=O)[OX2H1]", "smarts"),
    ("c[OH]", "smarts"),
    ("[#7]", "smarts"),
    ("[R2]", "smarts"),  # fused ring atoms: sets few pattern bits, screen keeps most
]


def oracle(query: str, query_type: str) -> list[int]:
    """CIDs whose structure contains the query, by brute force."""
    q = service.parse_query(query, query_type)  # type: ignore[arg-type]
    return [cid for cid, smi in ALL if Chem.MolFromSmiles(smi).HasSubstructMatch(q)]


def fresh_home(tmp_path, sdfs):
    paths = IndexPaths(tmp_path / "home")
    srcs = [SdfSource(name=chunk_name(p), local_path=p) for p in sdfs]
    build_chunks(paths, srcs, workers=1)
    return paths


def test_pattern_layout_and_manifest(substructure_home):
    paths = IndexPaths(substructure_home)
    manifest = paths.load_manifest()
    assert manifest is not None
    assert manifest.substructure is not None
    sub = manifest.substructure
    assert sub.fp_type == "RDKitPattern"
    assert sub.n_compounds == N_GOOD
    assert sub.n_failed == 0
    assert sub.chunks == INDEXED_CHUNKS
    assert sub.memory_bytes() == N_GOOD * (2048 // 64 + 2) * 8
    assert paths.pattern_fingerprints.exists()
    for name in INDEXED_CHUNKS:
        assert paths.pattern_h5(name).exists()
        assert paths.pattern_stats(name).exists()
    # Pattern stats live in a subdirectory and never pass for chunk stats.
    assert [c.name for c in paths.completed_chunks()] == list(CHUNK_NAMES)
    assert not list(paths.index_dir.glob("*.tmp"))


@pytest.mark.parametrize(("query", "query_type"), QUERIES)
def test_matches_brute_force(substructure_home, query, query_type):
    result = service.substructure_search(
        query, query_type=query_type, max_hits=1000, home=substructure_home
    )
    assert [h.cid for h in result.hits] == oracle(query, query_type)
    assert result.exhaustive is True
    assert result.n_checked == result.n_candidates >= len(result.hits)
    assert result.n_indexed == N_GOOD
    assert result.warning is None


def test_screen_prunes_and_hits_carry_metadata(substructure_home):
    result = service.substructure_search("c1ccccc1C(=O)O", home=substructure_home)
    assert [h.cid for h in result.hits] == [243, 338, 2244, 4133]
    assert result.n_candidates < N_GOOD
    aspirin = result.hits[2]
    assert aspirin.iupac_name == "2-acetyloxybenzoic acid"
    assert aspirin.pubchem_url.endswith("/2244")
    assert "similarity" not in aspirin.model_dump()


def test_max_hits_stops_early(substructure_home):
    everything = oracle("[#6]", "smarts")
    result = service.substructure_search(
        "[#6]", query_type="smarts", max_hits=3, home=substructure_home
    )
    assert [h.cid for h in result.hits] == everything[:3]
    assert result.n_checked == 3
    assert result.exhaustive is False


def test_max_candidates_caps_verification(substructure_home):
    result = service.substructure_search(
        "[#6]", query_type="smarts", max_candidates=2, home=substructure_home
    )
    assert len(result.hits) == 2
    assert result.n_checked == 2
    assert result.n_candidates == N_GOOD
    assert result.exhaustive is False


def test_disk_mode_matches_memory(substructure_home):
    query = service.parse_query("C(=O)O")
    mem = LocalIndex(substructure_home, mode="memory").substructure(query, 50, 1000)
    disk = LocalIndex(substructure_home, mode="disk").substructure(query, 50, 1000)
    assert [r["cid"] for r in disk.rows] == [r["cid"] for r in mem.rows]
    assert disk.n_candidates == mem.n_candidates


@pytest.mark.parametrize(
    ("query", "query_type", "error"),
    [
        ("C1CC", "smiles", MoleculeValidationError),
        ("", "smiles", MoleculeValidationError),
        ("[C", "smarts", MoleculeValidationError),
        ("  ", "smarts", MoleculeValidationError),
        ("CCO", "inchi", PubChemFinderError),
    ],
)
def test_invalid_queries(substructure_home, query, query_type, error):
    with pytest.raises(error):
        service.substructure_search(
            query, query_type=query_type, home=substructure_home
        )


def test_invalid_limits(substructure_home):
    with pytest.raises(PubChemFinderError, match="max_hits"):
        service.substructure_search("C", max_hits=0, home=substructure_home)
    with pytest.raises(PubChemFinderError, match="max_candidates"):
        service.substructure_search("C", max_candidates=0, home=substructure_home)


def test_index_without_substructure_explains_how_to_add_it(built_home):
    with pytest.raises(IndexNotBuilt, match="--with-substructure"):
        service.substructure_search("c1ccccc1", home=built_home)
    status = service.index_status(built_home)
    assert status.substructure is False
    assert status.substructure_size_bytes == 0
    assert "--with-substructure" in (status.hint or "")


def test_no_index_at_all(tmp_path):
    with pytest.raises(IndexNotBuilt, match="index build"):
        service.substructure_search("c1ccccc1", home=tmp_path)


def test_status_reports_substructure(substructure_home):
    status = service.index_status(substructure_home)
    assert status.substructure is True
    assert status.substructure_size_bytes > 0
    assert status.substructure_memory_bytes == N_GOOD * 34 * 8
    assert status.hint is None


def test_substructure_lifecycle(tmp_path, sample_sdfs):
    """Add to an existing index, keep it on for new chunks, then remove it."""
    paths = fresh_home(tmp_path, sample_sdfs[:1])
    finalize(paths)
    aspirin = "CC(=O)Oc1ccccc1C(=O)O"
    before = paths.fingerprints.read_bytes(), paths.fingerprints.stat().st_mtime_ns
    similar = service.similarity_search(aspirin, threshold=0.3, home=paths.home)

    # Adding substructure search leaves the similarity database alone.
    manifest = finalize(paths, substructure=True)
    assert manifest.substructure is not None
    assert manifest.substructure.chunks == [CHUNK_NAMES[0]]
    assert (
        paths.fingerprints.read_bytes(),
        paths.fingerprints.stat().st_mtime_ns,
    ) == before
    assert service.similarity_search(aspirin, threshold=0.3, home=paths.home) == similar

    # Nothing changed: a second finalize leaves the manifest alone.
    written = paths.manifest.stat().st_mtime_ns
    assert finalize(paths) == manifest
    assert paths.manifest.stat().st_mtime_ns == written

    # A later build without the flag keeps substructure search up to date.
    new_chunk = SdfSource(name=chunk_name(sample_sdfs[1]), local_path=sample_sdfs[1])
    build_chunks(paths, [new_chunk])
    manifest = finalize(paths)
    assert manifest.substructure is not None
    assert manifest.substructure.chunks == INDEXED_CHUNKS
    caffeine_core = service.substructure_search("Cn1cnc2c1cncn2", home=paths.home)
    assert [h.cid for h in caffeine_core.hits] == [2519]

    # --without-substructure removes the pattern database and its chunks.
    manifest = finalize(paths, substructure=False)
    assert manifest.substructure is None
    assert not paths.pattern_fingerprints.exists()
    assert not paths.pattern_chunks.exists()
    assert paths.load_manifest() == manifest
    assert service.similarity_search(aspirin, threshold=0.3, home=paths.home).hits
    with pytest.raises(IndexNotBuilt, match="--with-substructure"):
        service.substructure_search("CCO", home=paths.home)


def test_manifest_without_substructure_field_loads(substructure_home):
    data = json.loads(IndexPaths(substructure_home).manifest.read_text())
    del data["substructure"]
    assert Manifest.model_validate(data).substructure is None


def test_rdkit_version_change_warns_then_rebuilds(tmp_path, substructure_home):
    home = tmp_path / "home"
    shutil.copytree(substructure_home, home)
    paths = IndexPaths(home)
    manifest = paths.load_manifest()
    assert manifest is not None
    assert manifest.substructure is not None
    old = manifest.substructure.model_copy(update={"rdkit_version": "2000.01.1"})
    paths.manifest.write_text(
        manifest.model_copy(update={"substructure": old}).model_dump_json()
    )
    stale = paths.pattern_stats(INDEXED_CHUNKS[0])
    stats = json.loads(stale.read_text())
    stale.write_text(json.dumps({**stats, "rdkit_version": "2000.01.1"}))

    result = service.substructure_search("c1ccccc1", home=home)
    assert "Rebuild" in (result.warning or "")
    assert [h.cid for h in result.hits] == oracle("c1ccccc1", "smiles")

    events: list[str] = []
    rebuilt = finalize(paths, events.append)
    assert rebuilt.substructure is not None
    assert rebuilt.substructure.rdkit_version != "2000.01.1"
    assert any(INDEXED_CHUNKS[0] in e for e in events)
    assert not any(
        INDEXED_CHUNKS[1] in e and "pattern fingerprints" in e for e in events
    )
    assert service.substructure_search("c1ccccc1", home=home).warning is None


def test_parallel_pattern_build(tmp_path, sample_sdfs):
    paths = fresh_home(tmp_path, sample_sdfs[:2])
    built = build_pattern_chunks(paths, INDEXED_CHUNKS, workers=2)
    assert sorted(s.name for s in built) == INDEXED_CHUNKS
    assert sum(s.n_records for s in built) == N_GOOD
    assert build_pattern_chunks(paths, INDEXED_CHUNKS, workers=2) == []


def test_unparseable_smiles_are_counted(tmp_path, sample_sdfs):
    import pyarrow.parquet as pq

    paths = fresh_home(tmp_path, sample_sdfs[:1])
    meta = paths.chunk_parquet(CHUNK_NAMES[0])
    table = pq.read_table(meta)
    smiles = table.column("smiles").to_pylist()
    smiles[0] = "not-a-smiles"
    smiles[1] = None
    pq.write_table(table.set_column(1, "smiles", [smiles]), meta)
    [stats] = build_pattern_chunks(paths, [CHUNK_NAMES[0]])
    assert (stats.n_records, stats.n_failed) == (len(CHUNK_A) - 2, 2)


def test_pattern_build_needs_parseable_smiles(tmp_path, sample_sdfs):
    import pyarrow.parquet as pq

    paths = fresh_home(tmp_path, sample_sdfs[:1])
    manifest = finalize(paths)
    meta = paths.chunk_parquet(CHUNK_NAMES[0])
    table = pq.read_table(meta)
    pq.write_table(table.set_column(1, "smiles", [["x"] * table.num_rows]), meta)
    with pytest.raises(IndexBuildError, match="No stored SMILES"):
        finalize(paths, substructure=True)
    assert paths.load_manifest() == manifest  # the existing index is still usable


def test_small_metadata_batches_give_the_same_answer(substructure_home):
    query = service.parse_query("[#6]", "smarts")
    found = LocalIndex(substructure_home).substructure(
        query, max_hits=1000, max_candidates=1000, max_batch=2
    )
    assert [r["cid"] for r in found.rows] == oracle("[#6]", "smarts")
    assert found.exhaustive is True


def test_pattern_chunks_rebuilt_when_fpsim2_changes(tmp_path, substructure_home):
    """FPSim2 refuses to merge files written by different versions."""
    home = tmp_path / "home"
    shutil.copytree(substructure_home, home)
    paths = IndexPaths(home)
    stale = paths.pattern_stats(INDEXED_CHUNKS[0])
    stale.write_text(
        json.dumps({**json.loads(stale.read_text()), "fpsim2_version": "0.0.1"})
    )
    # A chunk from 2.1.0 records no FPSim2 version; its .h5 says it is current.
    legacy = paths.pattern_stats(INDEXED_CHUNKS[1])
    stats = json.loads(legacy.read_text())
    del stats["fpsim2_version"]
    legacy.write_text(json.dumps(stats))

    events: list[str] = []
    finalize(paths, events.append, force=True)
    rebuilt = [e for e in events if "pattern fingerprints for" in e]
    assert len(rebuilt) == 1
    assert INDEXED_CHUNKS[0] in rebuilt[0]
    assert [h.cid for h in service.substructure_search("c1ccccc1", home=home).hits] == (
        oracle("c1ccccc1", "smiles")
    )
