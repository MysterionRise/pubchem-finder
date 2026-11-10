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
import uuid
from collections.abc import Callable, Iterator, Sequence
from concurrent.futures import (
    FIRST_COMPLETED,
    Future,
    ProcessPoolExecutor,
    as_completed,
    wait,
)
from dataclasses import dataclass, field, replace
from datetime import UTC, datetime
from importlib.metadata import version
from pathlib import Path
from typing import Any, Literal, Protocol

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
from ..data.sdf import SdfReader, SdfRecord
from .layout import (
    FORMAT_VERSION,
    PATTERN_FP_TYPE,
    ChunkStats,
    IndexPaths,
    Manifest,
    PatternChunkStats,
    SourceFile,
    SubstructureInfo,
    chunk_key,
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


class RecordReader(Protocol):
    """Compounds to index, counting the records it had to skip (see SdfReader)."""

    n_ok: int
    n_failed: int

    def __iter__(self) -> Iterator[SdfRecord]: ...


def build_chunk(
    sdf_path: str,
    home: str,
    fp_params: dict[str, Any] | None = None,
    source: SourceFile | None = None,
    name: str | None = None,
) -> ChunkStats:
    """
    Fingerprint one SDF file into a chunk (see :func:`write_chunk`), named after
    the file unless ``name`` says otherwise.

    Runs in worker processes, so it takes and returns plain picklable values.
    """
    rdBase.DisableLog("rdApp.*")
    return write_chunk(
        IndexPaths(Path(home)),
        name or chunk_name(sdf_path),
        SdfReader(sdf_path),
        fp_params,
        source,
    )


def write_chunk(
    paths: IndexPaths,
    name: str,
    reader: RecordReader,
    fp_params: dict[str, Any] | None = None,
    source: SourceFile | None = None,
) -> ChunkStats:
    """
    Fingerprint ``reader``'s compounds into ``<name>.h5`` + ``<name>.parquet``
    in one pass and write ``<name>.json`` last to mark the chunk complete.

    A chunk that already exists stays as it is until the new files are complete,
    so a failed or refused rebuild loses nothing. The swap then takes the old
    ``<name>.json`` away first: if it is interrupted, the chunk reads as unbuilt
    rather than as the old build with new files.
    """
    params = dict(fp_params or DEFAULT_FP_PARAMS)
    paths.chunks.mkdir(parents=True, exist_ok=True)
    h5_tmp = paths.chunk_h5(name).with_suffix(".h5.tmp")
    pq_tmp = paths.chunk_parquet(name).with_suffix(".parquet.tmp")
    start = time.perf_counter()

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

    paths.chunk_stats(name).unlink(missing_ok=True)
    paths.pattern_stats(name).unlink(missing_ok=True)  # derived from the old build
    if reader.n_ok:
        os.replace(h5_tmp, paths.chunk_h5(name))
    else:  # FPSim2 cannot read an empty database
        h5_tmp.unlink(missing_ok=True)
        paths.chunk_h5(name).unlink(missing_ok=True)
    os.replace(pq_tmp, paths.chunk_parquet(name))
    paths.pattern_h5(name).unlink(missing_ok=True)
    stats = ChunkStats(
        name=name,
        n_records=reader.n_ok,
        n_failed=reader.n_failed,
        seconds=round(time.perf_counter() - start, 2),
        fp_type=FP_TYPE,
        fp_params=params,
        rdkit_version=rdkit.__version__,
        fpsim2_version=version("fpsim2"),
        source=source,
        build_id=uuid.uuid4().hex,
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
        with pq.ParquetFile(paths.chunk_parquet(name)) as meta:
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

    def source(self) -> SourceFile | None:
        """The local file's identity, or None for a PubChem chunk."""
        return SourceFile.of(self.local_path) if self.local_path is not None else None

    def label(self) -> str:
        return str(self.local_path) if self.local_path is not None else self.name

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
    # Chunks replaced under their old name. finalize needs no telling: it compares
    # each chunk's build_id with the one its manifest recorded.
    rebuilt: list[str] = field(default_factory=list)


def build_chunks(
    paths: IndexPaths,
    sources: Sequence[SdfSource],
    workers: int = 1,
    fp_params: dict[str, Any] | None = None,
    on_event: EventFn | None = None,
) -> BuildReport:
    """
    Build every source that has no completed chunk yet, and rebuild local files
    that changed since their chunk was built.

    Downloads happen in this process, one at a time; fingerprinting runs in up to
    ``workers`` processes. At most ``2 * workers`` SDF files are on disk at once,
    and downloaded files are deleted as soon as their chunk is complete. A chunk is
    named after its file, so two files with one name are refused rather than
    written over each other.

    Every check runs before anything on disk changes, and a chunk being rebuilt
    is replaced only once its new files are complete (see :func:`write_chunk`).
    """
    params = dict(fp_params or DEFAULT_FP_PARAMS)
    emit = on_event or (lambda msg: logger.info(msg))
    completed = paths.completed_chunks()
    sources = _check_unique_names(sources, completed)
    _check_params(completed, paths, params)
    by_name = {stats.name: stats for stats in completed}
    report = BuildReport()
    todo: list[SdfSource] = []
    adopt: list[tuple[ChunkStats, SourceFile, str]] = []
    reasons: dict[str, str] = {}
    for src in sources:
        stats = by_name.get(src.name)
        if stats is None:
            todo.append(src)
            continue
        action, note = _plan(src, stats, paths.chunk_stats(src.name))
        if action == "rebuild":
            reasons[src.name] = note
            todo.append(src)
            continue
        report.skipped.append(src.name)
        current = src.source()
        if action == "adopt" and current is not None:
            adopt.append((stats, current, note))
    if todo:
        # Before any download or rebuild: FPSim2 can't merge new chunks with ones
        # another version wrote. Chunks about to be replaced don't count.
        _check_versions(paths, [c for c in completed if c.name not in reasons])
    paths.chunks.mkdir(parents=True, exist_ok=True)
    for stats, current, note in adopt:
        emit(note)
        write_json_atomic(
            paths.chunk_stats(stats.name), stats.model_copy(update={"source": current})
        )
    for name, note in reasons.items():
        emit(f"rebuilding {name}: {note}")
        report.rebuilt.append(name)
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
            try:
                stats = build_chunk(
                    str(path), str(paths.home), params, src.source(), src.name
                )
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
            fut = pool.submit(
                build_chunk, str(path), str(paths.home), params, src.source(), src.name
            )
            pending[fut] = (src, path)
            drain(2 * workers - 1)
        drain(0)
    return report


def _same_file(a: Path | str, b: Path | str) -> bool:
    try:
        return os.path.samefile(a, b)
    except OSError:  # one of them is gone
        return False


def _check_unique_names(
    sources: Sequence[SdfSource], completed: Sequence[ChunkStats]
) -> list[SdfSource]:
    """Refuse sources that would share a chunk, comparing names without case
    (``Library`` and ``library`` are the same files on macOS and Windows).

    Returns the sources, each under the spelling its existing chunk already has.
    """
    existing = {chunk_key(stats.name): stats for stats in completed}
    seen: dict[str, SdfSource] = {}
    out: list[SdfSource] = []
    for src in sources:
        key = chunk_key(src.name)
        other = seen.setdefault(key, src)
        if other is not src:
            if (
                other.local_path is not None
                and src.local_path is not None
                and _same_file(other.local_path, src.local_path)
            ):
                raise IndexBuildError(f"{src.label()} is given twice.")
            raise IndexBuildError(
                f"{other.label()} and {src.label()} would both be chunk {src.name!r} "
                "(names are compared without case): each file becomes a chunk named "
                "after it. Rename one of them."
            )
        have = existing.get(key)
        if have is not None and have.name != src.name:
            same = (
                have.source is not None
                and src.local_path is not None
                and _same_file(have.source.path, src.local_path)
            )
            if not same:
                raise IndexBuildError(
                    f"{src.label()} would be chunk {src.name!r}, and this index "
                    f"already has chunk {have.name!r}: names are compared without "
                    "case. Rename the file."
                )
            src = replace(src, name=have.name)  # the same file, typed in another case
        out.append(src)
    return out


def _plan(
    src: SdfSource, stats: ChunkStats, stats_path: Path
) -> tuple[Literal["skip", "adopt", "rebuild"], str]:
    """What to do about a source whose chunk is already complete, and why.

    ``adopt`` keeps the chunk and records ``src`` as its file: the same file at a
    new place, or the file behind a chunk an earlier version built.
    """
    current = src.source()
    if current is None:  # a PubChem chunk: its name pins its content
        return "skip", ""
    recorded = stats.source
    if recorded is None:  # 2.2.0 and earlier did not record the file
        if current.mtime_ns <= stats_path.stat().st_mtime_ns:
            return "adopt", f"{src.name}: recorded {current.path} as its file"
        return "rebuild", f"{current.path} is newer than its chunk"
    unchanged = (recorded.size, recorded.mtime_ns) == (current.size, current.mtime_ns)
    if recorded.path == current.path or _same_file(recorded.path, current.path):
        if unchanged:
            return "skip", ""
        return "rebuild", f"{current.path} changed since it was indexed"
    if unchanged:
        return "adopt", f"{src.name}: its file moved to {current.path}"
    raise IndexBuildError(
        f"Chunk {src.name!r} already holds {recorded.path}, and {current.path} is a "
        "different file that would replace it: each file becomes a chunk named after "
        "it. Rename one of them."
    )


def _check_params(
    completed: Sequence[ChunkStats], paths: IndexPaths, params: dict[str, Any]
) -> None:
    for stats in completed:
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

    If the existing index already holds exactly the completed chunks, in the builds
    now on disk, it is left alone (merging all of PubChem takes a while) unless
    ``force`` is set. That is read from disk, not from what this process did: a
    chunk rebuilt by a run that then failed is merged by the next one.
    ``substructure`` turns the pattern-fingerprint database on or off; ``None``
    keeps whatever the index has, or was being given when a merge was interrupted
    (its pattern chunks are the record of that). Turning it on builds any missing
    pattern chunks (in up to ``workers`` processes) without touching the
    similarity database; only turning it off removes pattern data.
    """
    emit = on_event or (lambda msg: logger.info(msg))
    chunks = [c for c in paths.completed_chunks() if c.n_records > 0]
    if not chunks:
        raise IndexBuildError(f"No built chunks in {paths.chunks}; nothing to finalize")
    names = [c.name for c in chunks]
    current = paths.load_manifest()
    merged = current is not None and _holds(current, chunks)
    old_sub = current.substructure if current is not None else None
    if substructure is None:
        # The manifest goes without substructure search while the pattern database
        # is merged, so after an interrupted merge it is the pattern chunks that
        # say the user asked for it.
        want_sub = old_sub is not None or bool(paths.completed_pattern_chunks())
    else:
        want_sub = substructure
    sub_fresh = (
        not force
        and merged
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

    if not force and current is not None and merged:
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
        if substructure is False and _remove_pattern_data(paths):
            emit("substructure search removed")
        emit(f"index up to date: {manifest.n_compounds:,} compounds")
        return manifest
    manifest = manifest.model_copy(update={"substructure": sub})
    write_json_atomic(paths.manifest, manifest)
    if substructure is False and _remove_pattern_data(paths):
        emit("substructure search removed")
    extra = " with substructure search" if sub is not None else ""
    emit(f"index ready: {manifest.n_compounds:,} compounds{extra}")
    return manifest


def _remove_pattern_data(paths: IndexPaths) -> bool:
    """Delete the pattern database and its chunks; whether there was anything."""
    found = paths.pattern_fingerprints.exists() or paths.pattern_chunks.is_dir()
    paths.pattern_fingerprints.unlink(missing_ok=True)
    shutil.rmtree(paths.pattern_chunks, ignore_errors=True)
    return found


def _holds(manifest: Manifest, chunks: Sequence[ChunkStats]) -> bool:
    """Whether the merged index holds exactly ``chunks``, in the builds on disk."""
    if manifest.chunks != [c.name for c in chunks]:
        return False
    if manifest.chunk_builds is None:
        # Merged by 2.2.1 or earlier, which kept no record of builds. 2.2.1 could
        # rebuild an --sdf chunk and fail before merging it, and such chunks are
        # the ones with a recorded file, so an index with one is merged once more.
        return all(c.source is None for c in chunks)
    return manifest.chunk_builds == {c.name: c.build_id for c in chunks}


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
        format_version=FORMAT_VERSION,
        fp_type=chunks[0].fp_type,
        fp_params=params,
        n_compounds=sum(c.n_records for c in chunks),
        n_failed=sum(c.n_failed for c in all_stats),
        chunks=[c.name for c in chunks],
        built_at=datetime.now(UTC),
        pubchem_finder_version=__version__,
        rdkit_version=rdkit.__version__,
        fpsim2_version=version("fpsim2"),
        chunk_builds={c.name: c.build_id for c in chunks},
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
    _merge_h5(
        files,
        paths.pattern_fingerprints,
        on_failure="Similarity search is unaffected; run the build again to finish "
        "substructure search.",
    )
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
    on_failure: str = "The current index is unchanged.",
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
            f"Merging into {dest} failed: {exc}. {on_failure}"
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
