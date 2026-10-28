"""Build the local similarity index from PubChem SDF files.

Each SDF file becomes one *chunk*: an FPSim2 fingerprint file and a Parquet
metadata file, written in a single pass over the SDF. Chunks are independent,
so they build in parallel and an interrupted build resumes where it stopped.
Substructure search is opt-in: :func:`build_pattern_chunk` derives RDKit
pattern fingerprints from a chunk's Parquet SMILES (no SDF needed), so the
screen and the verification step see the same molecules.
:func:`finalize` merges all completed chunks into one popcount-sorted FPSim2
database (and, if enabled, a second one of pattern fingerprints) and writes the
manifest that marks the index as ready.
"""

import logging
import multiprocessing
import os
import shutil
import time
from collections.abc import Callable, Iterator, Sequence
from concurrent.futures import (
    FIRST_COMPLETED,
    Future,
    ProcessPoolExecutor,
    as_completed,
    wait,
)
from dataclasses import dataclass, field
from datetime import UTC, datetime
from importlib.metadata import version
from pathlib import Path
from typing import Any

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import rdkit
import tables as tb
from FPSim2.io import create_db_file
from FPSim2.io.backends.pytables import merge_db_files
from rdkit import Chem, rdBase

from .. import __version__
from ..core.exceptions import IndexBuildError, PubChemFinderError
from ..data.sdf import SdfReader
from .layout import (
    PATTERN_FP_TYPE,
    ChunkStats,
    IndexPaths,
    Manifest,
    PatternChunkStats,
    SubstructureInfo,
    chunk_name,
    write_json_atomic,
)

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
POPCNT_BLOCK = 1 << 22  # rows of popcount read at a time when sorting

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
        rdkit_version=rdkit.__version__,
        fpsim2_version=version("fpsim2"),
    )
    write_json_atomic(paths.chunk_stats(name), stats)
    return stats


def build_pattern_chunk(home: str, name: str) -> PatternChunkStats:
    """
    Pattern-fingerprint one built chunk from the SMILES in ``<chunk>.parquet``
    into ``pattern/<chunk>.h5``, writing ``pattern/<chunk>.json`` last.

    Runs in worker processes, so it takes and returns plain picklable values.
    """
    rdBase.DisableLog("rdApp.*")
    paths = IndexPaths(Path(home))
    paths.pattern_chunks.mkdir(parents=True, exist_ok=True)
    h5_tmp = paths.pattern_h5(name).with_suffix(".h5.tmp")
    start = time.perf_counter()
    counts = {"ok": 0, "failed": 0}

    def mols() -> Iterator[tuple[Chem.Mol, int]]:
        meta = pq.ParquetFile(paths.chunk_parquet(name))
        for batch in meta.iter_batches(
            batch_size=METADATA_BATCH, columns=["cid", "smiles"]
        ):
            cids = batch.column("cid").to_pylist()
            smiles = batch.column("smiles").to_pylist()
            for cid, smi in zip(cids, smiles, strict=True):
                mol = Chem.MolFromSmiles(smi) if smi else None
                if mol is None:
                    counts["failed"] += 1
                    continue
                counts["ok"] += 1
                yield mol, cid

    records = _peek(mols())
    if records is not None:
        create_db_file(
            records, str(h5_tmp), "rdkit", PATTERN_FP_TYPE, {}, sort_by_popcnt=False
        )
        os.replace(h5_tmp, paths.pattern_h5(name))
    stats = PatternChunkStats(
        name=name,
        n_records=counts["ok"],
        n_failed=counts["failed"],
        seconds=round(time.perf_counter() - start, 2),
        rdkit_version=rdkit.__version__,
        fpsim2_version=version("fpsim2"),
    )
    write_json_atomic(paths.pattern_stats(name), stats)
    return stats


def build_pattern_chunks(
    paths: IndexPaths,
    names: Sequence[str],
    workers: int = 1,
    on_event: EventFn | None = None,
) -> list[PatternChunkStats]:
    """Build pattern fingerprints for every chunk in ``names`` that lacks them."""
    emit = on_event or (lambda msg: logger.info(msg))
    done = paths.completed_pattern_chunks()
    # Pattern fingerprints can change between RDKit releases, and a screen built
    # with different bits than the query would miss matches; FPSim2 also refuses to
    # merge files written by different versions: rebuild stale ones.
    todo = [n for n in names if n not in done or _pattern_stale(paths, done[n])]
    built: list[PatternChunkStats] = []

    def report(stats: PatternChunkStats) -> None:
        built.append(stats)
        emit(
            f"pattern fingerprints for {stats.name}: {stats.n_records:,} compounds"
            f" ({stats.n_failed:,} skipped) in {stats.seconds:.0f}s"
        )

    if workers <= 1 or len(todo) <= 1:
        for name in todo:
            report(build_pattern_chunk(str(paths.home), name))
        return built
    ctx = multiprocessing.get_context("spawn")
    with ProcessPoolExecutor(max_workers=workers, mp_context=ctx) as pool:
        futures = {
            pool.submit(build_pattern_chunk, str(paths.home), name): name
            for name in todo
        }
        for fut in as_completed(futures):
            try:
                report(fut.result())
            except Exception as exc:
                raise IndexBuildError(
                    f"Pattern fingerprints for {futures[fut]} failed: {exc}"
                ) from exc
    return built


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
    if todo:
        # Before any download: FPSim2 can't merge new chunks with old ones.
        _check_versions(paths, paths.completed_chunks())

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
            try:
                stats = build_chunk(str(path), str(paths.home), params)
            except PubChemFinderError:
                raise
            except Exception as exc:
                raise IndexBuildError(f"Building {src.name} failed: {exc}") from exc
            done(src, path, stats)
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


def runtime_versions() -> tuple[str, str]:
    """The RDKit and FPSim2 versions this process fingerprints with."""
    return rdkit.__version__, version("fpsim2")


def h5_versions(path: Path) -> tuple[str, str]:
    """The RDKit and FPSim2 versions recorded in an FPSim2 file, which FPSim2
    compares (with the fingerprint settings) before it merges files."""
    with tb.open_file(str(path), mode="r") as f:
        return str(f.root.config[2]), str(f.root.config[3])


def _chunk_versions(paths: IndexPaths, stats: ChunkStats) -> tuple[str, str] | None:
    if stats.n_records == 0:
        return None  # no .h5, nothing to merge
    if stats.rdkit_version and stats.fpsim2_version:
        return stats.rdkit_version, stats.fpsim2_version
    return h5_versions(paths.chunk_h5(stats.name))


def _pattern_stale(paths: IndexPaths, stats: PatternChunkStats) -> bool:
    rdkit_now, fpsim2_now = runtime_versions()
    if stats.rdkit_version != rdkit_now:
        return True
    if stats.n_records == 0:
        return False
    built_with = stats.fpsim2_version or h5_versions(paths.pattern_h5(stats.name))[1]
    return built_with != fpsim2_now


def _version_help(paths: IndexPaths) -> str:
    return (
        "FPSim2 can't merge fingerprints written by different versions. The current "
        "index is unchanged and still searchable. To grow it, build the larger index "
        "in a new data directory (`pubchem-finder index build --limit N --home NEW_DIR`)"
        f" and point PUBCHEM_FINDER_HOME at it once it's ready, or delete {paths.home}"
        " and build again."
    )


def _check_versions(paths: IndexPaths, chunks: Sequence[ChunkStats]) -> None:
    """Refuse to add chunks next to ones fingerprinted by another RDKit/FPSim2."""
    now = runtime_versions()
    old: dict[tuple[str, str], int] = {}
    for stats in chunks:
        built_with = _chunk_versions(paths, stats)
        if built_with is not None and built_with != now:
            old[built_with] = old.get(built_with, 0) + 1
    if old:
        rdkit_then, fpsim2_then = next(iter(old))
        raise IndexBuildError(
            f"{sum(old.values())} chunk(s) in {paths.chunks} were fingerprinted with "
            f"RDKit {rdkit_then} / FPSim2 {fpsim2_then}; this pubchem-finder runs RDKit "
            f"{now[0]} / FPSim2 {now[1]}. {_version_help(paths)}"
        )


def finalize(
    paths: IndexPaths,
    on_event: EventFn | None = None,
    force: bool = False,
    substructure: bool | None = None,
    workers: int = 1,
) -> Manifest:
    """Merge all completed chunks into the searchable index and write its manifest.

    If the existing index already contains exactly the completed chunks, it is
    left alone (merging all of PubChem takes a while) unless ``force`` is set.
    ``substructure`` turns the pattern-fingerprint database on or off; ``None``
    keeps whatever the current index has. Turning it on builds any missing
    pattern chunks (in up to ``workers`` processes) without touching the
    similarity database.
    """
    emit = on_event or (lambda msg: logger.info(msg))
    chunks = [c for c in paths.completed_chunks() if c.n_records > 0]
    if not chunks:
        raise IndexBuildError(f"No built chunks in {paths.chunks}; nothing to finalize")
    names = [c.name for c in chunks]
    current = paths.load_manifest()
    old_sub = current.substructure if current is not None else None
    want_sub = old_sub is not None if substructure is None else substructure
    sub_fresh = (
        not force
        and old_sub is not None
        and old_sub.chunks == names
        and old_sub.rdkit_version == rdkit.__version__
        and paths.pattern_fingerprints.exists()
    )
    pattern_files: list[str] = []
    if want_sub and not sub_fresh:
        # Before touching either database, so a failure leaves the current
        # index usable.
        build_pattern_chunks(paths, names, workers, emit)
        pattern_files = _pattern_files(paths, names)

    if not force and current is not None and current.chunks == names:
        manifest = current
        similarity_rebuilt = False
    else:
        manifest = _merge_similarity(paths, chunks, emit)
        similarity_rebuilt = True
        # Similarity search works again while the pattern database is merged; the
        # manifest gains substructure search below, once that merge is done.
        write_json_atomic(paths.manifest, manifest)

    sub: SubstructureInfo | None = None
    if want_sub and sub_fresh:
        sub = old_sub
    elif want_sub:
        if not similarity_rebuilt and old_sub is not None:
            # Running servers must not pair the old popcount bins with the new
            # file: withdraw substructure search until the merge is done.
            write_json_atomic(
                paths.manifest, manifest.model_copy(update={"substructure": None})
            )
        sub = _merge_patterns(paths, names, pattern_files, emit)

    if not similarity_rebuilt and sub == old_sub:
        emit(f"index up to date: {manifest.n_compounds:,} compounds")
        return manifest
    manifest = manifest.model_copy(update={"substructure": sub})
    write_json_atomic(paths.manifest, manifest)
    if sub is None and paths.pattern_fingerprints.exists():
        paths.pattern_fingerprints.unlink()
        shutil.rmtree(paths.pattern_chunks, ignore_errors=True)
        emit("substructure search removed")
    extra = " with substructure search" if sub is not None else ""
    emit(f"index ready: {manifest.n_compounds:,} compounds{extra}")
    return manifest


def _merge_similarity(
    paths: IndexPaths, chunks: list[ChunkStats], emit: EventFn
) -> Manifest:
    params = chunks[0].fp_params
    if any(c.fp_params != params or c.fp_type != chunks[0].fp_type for c in chunks):
        raise IndexBuildError("Chunks were built with different fingerprint settings")
    paths.index_dir.mkdir(parents=True, exist_ok=True)
    files = [str(paths.chunk_h5(c.name)) for c in chunks]
    emit(f"merging {len(files)} chunk(s) into {paths.fingerprints}")
    # The old manifest stays until the merged file is complete, so a failed merge
    # leaves the previous index searchable; it is withdrawn only for the swap.
    _merge_h5(files, paths.fingerprints, withdraw=paths.manifest, home=paths.home)
    all_stats = paths.completed_chunks()
    return Manifest(
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


def _pattern_files(paths: IndexPaths, names: list[str]) -> list[str]:
    done = paths.completed_pattern_chunks()
    missing = [n for n in names if n not in done]
    if missing:
        raise IndexBuildError(f"Pattern fingerprints missing for {', '.join(missing)}")
    files = [str(paths.pattern_h5(n)) for n in names if done[n].n_records > 0]
    if not files:
        raise IndexBuildError(
            "No stored SMILES could be parsed for pattern fingerprints"
        )
    return files


def _merge_patterns(
    paths: IndexPaths, names: list[str], files: list[str], emit: EventFn
) -> SubstructureInfo:
    stats = [paths.completed_pattern_chunks()[n] for n in names]
    emit(f"merging {len(files)} pattern chunk(s) into {paths.pattern_fingerprints}")
    _merge_h5(files, paths.pattern_fingerprints)
    return SubstructureInfo(
        n_compounds=sum(s.n_records for s in stats),
        n_failed=sum(s.n_failed for s in stats),
        chunks=list(names),
        built_at=datetime.now(UTC),
        rdkit_version=rdkit.__version__,
    )


def _merge_h5(
    files: list[str],
    dest: Path,
    withdraw: Path | None = None,
    home: Path | None = None,
) -> None:
    """Merge FPSim2 files into one popcount-sorted file at ``dest``, atomically.

    Nothing at ``dest`` changes unless the merged file is complete. ``withdraw``
    (the manifest) is moved aside while ``dest`` is replaced, so running servers
    never pair it with a half-written file, and put back if the swap fails.
    """
    dest.parent.mkdir(parents=True, exist_ok=True)
    unsorted = dest.with_name(dest.name + ".unsorted")
    tmp = dest.with_name(dest.name + ".tmp")
    for path in (unsorted, tmp):
        path.unlink(missing_ok=True)
    _check_mergeable(files, home or dest.parent)
    try:
        if len(files) == 1:
            shutil.copyfile(files[0], unsorted)
        else:
            merge_db_files(files, str(unsorted), sort_by_popcnt=False)
        sort_fpsim2_file(unsorted, tmp)
    except Exception as exc:
        tmp.unlink(missing_ok=True)
        if isinstance(exc, PubChemFinderError):
            raise
        raise IndexBuildError(
            f"Merging into {dest} failed: {exc}. The current index is unchanged."
        ) from exc
    finally:
        unsorted.unlink(missing_ok=True)
    _swap_in(tmp, dest, withdraw)


def _check_mergeable(files: list[str], home: Path) -> None:
    """FPSim2's merge raises a bare ValueError on mixed versions; say why first."""
    if len(files) < 2:
        return
    found: dict[tuple[str, str], list[str]] = {}
    for file in files:
        found.setdefault(h5_versions(Path(file)), []).append(Path(file).stem)
    if len(found) > 1:
        groups = "; ".join(
            f"RDKit {r} / FPSim2 {f}: {', '.join(names[:3])}"
            + (f" and {len(names) - 3} more" if len(names) > 3 else "")
            for (r, f), names in found.items()
        )
        raise IndexBuildError(
            f"Chunks were fingerprinted by different versions ({groups}). "
            f"{_version_help(IndexPaths(home))}"
        )


def _swap_in(tmp: Path, dest: Path, withdraw: Path | None) -> None:
    aside = withdraw.with_name(withdraw.name + ".old") if withdraw else None
    if withdraw is not None and aside is not None and withdraw.exists():
        os.replace(withdraw, aside)
    else:
        aside = None
    try:
        os.replace(tmp, dest)
    except OSError as exc:
        if withdraw is not None and aside is not None:
            os.replace(aside, withdraw)
        tmp.unlink(missing_ok=True)
        raise IndexBuildError(
            f"Could not replace {dest}: {exc}. If a server is using the index, stop "
            "it and run the build again; the current index is unchanged."
        ) from exc
    if aside is not None:
        aside.unlink(missing_ok=True)


def sort_fpsim2_file(src: Path, dest: Path) -> None:
    """
    Write a popcount-sorted copy of the FPSim2 database ``src`` to ``dest``.

    Produces the same file as FPSim2's ``sort_db_file``, whose popcount-bin step
    iterates every row in Python: ~30 us per row for 1024-bit Morgan, ~2 ms per
    row for 2048-bit pattern fingerprints (days for all of PubChem). Here the
    bins come from the sorted popcount column with numpy, in blocks.
    """
    filters = tb.Filters(complib="blosc2", complevel=9, fletcher32=False)
    with (
        tb.open_file(str(src), mode="r") as source,
        tb.open_file(str(dest), mode="w") as target,
    ):
        fp_type, fp_params = source.root.config[0], source.root.config[1]
        fps = source.root.fps.copy(
            target.root,
            "fps",
            filters=filters,
            copyuserattrs=True,
            overwrite=True,
            chunkshape="auto",
            sortby="popcnt",
            check_CSI=True,
            propindexes=True,
        )
        bins: dict[int, list[int]] = {}
        for start in range(0, fps.nrows, POPCNT_BLOCK):
            block = fps.read(start, start + POPCNT_BLOCK, field="popcnt")
            values, first, counts = np.unique(
                block, return_index=True, return_counts=True
            )
            for value, offset, count in zip(values, first, counts, strict=True):
                lo_hi = bins.setdefault(int(value), [start + int(offset), 0])
                lo_hi[1] = start + int(offset) + int(count)
        config = target.create_vlarray(target.root, "config", atom=tb.ObjectAtom())
        config.append(fp_type)
        config.append(fp_params)
        config.append(rdkit.__version__)
        config.append(version("fpsim2"))
        config.append([(value, (lo, hi)) for value, (lo, hi) in sorted(bins.items())])
