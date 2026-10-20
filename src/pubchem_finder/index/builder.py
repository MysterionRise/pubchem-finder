"""Build the local similarity index from PubChem SDF files.

Each SDF file becomes one *chunk*: an FPSim2 fingerprint file and a Parquet
metadata file, written in a single pass over the SDF. Chunks are independent,
so they build in parallel and an interrupted build resumes where it stopped.
:func:`finalize` merges all completed chunks into one popcount-sorted FPSim2
database and writes the manifest that marks the index as ready.
"""

import logging
import multiprocessing
import os
import shutil
import time
from collections.abc import Callable, Iterator, Sequence
from concurrent.futures import FIRST_COMPLETED, Future, ProcessPoolExecutor, wait
from dataclasses import dataclass, field
from datetime import UTC, datetime
from importlib.metadata import version
from pathlib import Path
from typing import Any

import pyarrow as pa
import pyarrow.parquet as pq
import rdkit
from FPSim2.io import create_db_file
from FPSim2.io.backends.pytables import merge_db_files, sort_db_file
from rdkit import Chem, rdBase

from .. import __version__
from ..core.exceptions import IndexBuildError
from ..data.sdf import SdfReader
from .layout import ChunkStats, IndexPaths, Manifest, chunk_name, write_json_atomic

logger = logging.getLogger(__name__)

FP_TYPE = "Morgan"
DEFAULT_FP_PARAMS: dict[str, Any] = {"radius": 2, "fpSize": 1024}

METADATA_SCHEMA = pa.schema(
    [
        ("cid", pa.int64()),
        ("smiles", pa.string()),
        ("inchikey", pa.string()),
        ("iupac_name", pa.string()),
        ("molecular_formula", pa.string()),
        ("molecular_weight", pa.float64()),
    ]
)
METADATA_BATCH = 50_000

EventFn = Callable[[str], None]


def build_chunk(
    sdf_path: str, home: str, fp_params: dict[str, Any] | None = None
) -> ChunkStats:
    """
    Fingerprint one SDF file into ``<chunk>.h5`` + ``<chunk>.parquet`` and
    write ``<chunk>.json`` last to mark it complete.

    Runs in worker processes, so it takes and returns plain picklable values.
    """
    rdBase.DisableLog("rdApp.*")
    params = dict(fp_params or DEFAULT_FP_PARAMS)
    paths = IndexPaths(Path(home))
    paths.chunks.mkdir(parents=True, exist_ok=True)
    name = chunk_name(sdf_path)
    h5_tmp = paths.chunk_h5(name).with_suffix(".h5.tmp")
    pq_tmp = paths.chunk_parquet(name).with_suffix(".parquet.tmp")
    start = time.perf_counter()
    reader = SdfReader(sdf_path)

    with pq.ParquetWriter(pq_tmp, METADATA_SCHEMA, compression="zstd") as writer:
        rows: dict[str, list[Any]] = {col: [] for col in METADATA_SCHEMA.names}

        def flush() -> None:
            if rows["cid"]:
                writer.write_table(pa.table(rows, schema=METADATA_SCHEMA))
                for col in rows.values():
                    col.clear()

        def mols() -> Iterator[tuple[Chem.Mol, int]]:
            for rec in reader:
                rows["cid"].append(rec.cid)
                rows["smiles"].append(rec.smiles)
                rows["inchikey"].append(rec.inchikey)
                rows["iupac_name"].append(rec.iupac_name)
                rows["molecular_formula"].append(rec.molecular_formula)
                rows["molecular_weight"].append(rec.molecular_weight)
                if len(rows["cid"]) >= METADATA_BATCH:
                    flush()
                yield rec.mol, rec.cid

        records = _peek(mols())
        if records is not None:
            # FPSim2 mutates fp_params; give it a copy.
            create_db_file(
                records,
                str(h5_tmp),
                "rdkit",
                FP_TYPE,
                dict(params),
                sort_by_popcnt=False,
            )
        flush()

    if reader.n_ok:
        os.replace(h5_tmp, paths.chunk_h5(name))
    else:
        h5_tmp.unlink(missing_ok=True)  # FPSim2 cannot read an empty database
    os.replace(pq_tmp, paths.chunk_parquet(name))
    stats = ChunkStats(
        name=name,
        n_records=reader.n_ok,
        n_failed=reader.n_failed,
        seconds=round(time.perf_counter() - start, 2),
        fp_type=FP_TYPE,
        fp_params=params,
    )
    write_json_atomic(paths.chunk_stats(name), stats)
    return stats


def _peek(it: Iterator[tuple[Chem.Mol, int]]) -> Iterator[tuple[Chem.Mol, int]] | None:
    """Return an equivalent iterator, or None if ``it`` is empty."""
    try:
        first = next(it)
    except StopIteration:
        return None

    def chained() -> Iterator[tuple[Chem.Mol, int]]:
        yield first
        yield from it

    return chained()


@dataclass
class SdfSource:
    """An SDF file to index: either local, or a PubChem chunk fetched on demand."""

    name: str
    local_path: Path | None = None
    fetch: Callable[[], Path] | None = None
    delete_after: bool = False

    def path(self) -> Path:
        if self.local_path is not None:
            return self.local_path
        if self.fetch is None:
            raise IndexBuildError(f"No way to obtain {self.name}")
        return self.fetch()


@dataclass
class BuildReport:
    built: list[ChunkStats] = field(default_factory=list)
    skipped: list[str] = field(default_factory=list)


def build_chunks(
    paths: IndexPaths,
    sources: Sequence[SdfSource],
    workers: int = 1,
    fp_params: dict[str, Any] | None = None,
    on_event: EventFn | None = None,
) -> BuildReport:
    """
    Build every source that has no completed chunk yet.

    Downloads happen in this process, one at a time; fingerprinting runs in up to
    ``workers`` processes. At most ``2 * workers`` SDF files are on disk at once,
    and downloaded files are deleted as soon as their chunk is complete.
    """
    params = dict(fp_params or DEFAULT_FP_PARAMS)
    emit = on_event or (lambda msg: logger.info(msg))
    paths.chunks.mkdir(parents=True, exist_ok=True)
    _check_params(paths, params)
    report = BuildReport()
    todo = []
    for src in sources:
        if paths.chunk_stats(src.name).exists():
            report.skipped.append(src.name)
        else:
            todo.append(src)
    if report.skipped:
        emit(f"{len(report.skipped)} chunk(s) already built, skipping")

    def done(src: SdfSource, path: Path, stats: ChunkStats) -> None:
        report.built.append(stats)
        emit(
            f"built {stats.name}: {stats.n_records:,} compounds"
            f" ({stats.n_failed:,} skipped) in {stats.seconds:.0f}s"
        )
        if src.delete_after:
            path.unlink(missing_ok=True)

    if workers <= 1:
        for src in todo:
            path = src.path()
            done(src, path, build_chunk(str(path), str(paths.home), params))
        return report

    ctx = multiprocessing.get_context("spawn")
    with ProcessPoolExecutor(max_workers=workers, mp_context=ctx) as pool:
        pending: dict[Future[ChunkStats], tuple[SdfSource, Path]] = {}

        def drain(limit: int) -> None:
            while len(pending) > limit:
                finished, _ = wait(pending, return_when=FIRST_COMPLETED)
                for fut in finished:
                    src, path = pending.pop(fut)
                    try:
                        done(src, path, fut.result())
                    except Exception as exc:
                        raise IndexBuildError(
                            f"Building {src.name} failed: {exc}"
                        ) from exc

        for src in todo:
            path = src.path()
            emit(f"fingerprinting {src.name}")
            fut = pool.submit(build_chunk, str(path), str(paths.home), params)
            pending[fut] = (src, path)
            drain(2 * workers - 1)
        drain(0)
    return report


def _check_params(paths: IndexPaths, params: dict[str, Any]) -> None:
    for stats in paths.completed_chunks():
        if stats.fp_type != FP_TYPE or stats.fp_params != params:
            raise IndexBuildError(
                f"Existing chunks in {paths.chunks} use {stats.fp_type} {stats.fp_params}, "
                f"not {FP_TYPE} {params}. Use a different --home or delete the chunks."
            )


def finalize(
    paths: IndexPaths, on_event: EventFn | None = None, force: bool = False
) -> Manifest:
    """Merge all completed chunks into the searchable index and write its manifest.

    If the existing index already contains exactly the completed chunks, it is
    left alone (merging all of PubChem takes a while) unless ``force`` is set.
    """
    emit = on_event or (lambda msg: logger.info(msg))
    chunks = [c for c in paths.completed_chunks() if c.n_records > 0]
    if not chunks:
        raise IndexBuildError(f"No built chunks in {paths.chunks}; nothing to finalize")
    current = paths.load_manifest()
    if not force and current is not None and current.chunks == [c.name for c in chunks]:
        emit(f"index up to date: {current.n_compounds:,} compounds")
        return current
    params = chunks[0].fp_params
    if any(c.fp_params != params or c.fp_type != chunks[0].fp_type for c in chunks):
        raise IndexBuildError("Chunks were built with different fingerprint settings")

    paths.index_dir.mkdir(parents=True, exist_ok=True)
    paths.manifest.unlink(missing_ok=True)  # the index is not ready until we finish
    tmp = paths.fingerprints.with_suffix(".h5.tmp")
    tmp.unlink(missing_ok=True)
    files = [str(paths.chunk_h5(c.name)) for c in chunks]
    emit(f"merging {len(files)} chunk(s) into {paths.fingerprints}")
    if len(files) == 1:
        shutil.copyfile(files[0], tmp)
        sort_db_file(str(tmp))
    else:
        merge_db_files(files, str(tmp), sort_by_popcnt=True)
    os.replace(tmp, paths.fingerprints)

    all_stats = paths.completed_chunks()
    manifest = Manifest(
        fp_type=chunks[0].fp_type,
        fp_params=params,
        n_compounds=sum(c.n_records for c in chunks),
        n_failed=sum(c.n_failed for c in all_stats),
        chunks=[c.name for c in chunks],
        built_at=datetime.now(UTC),
        pubchem_finder_version=__version__,
        rdkit_version=rdkit.__version__,
        fpsim2_version=version("fpsim2"),
    )
    write_json_atomic(paths.manifest, manifest)
    emit(f"index ready: {manifest.n_compounds:,} compounds")
    return manifest
