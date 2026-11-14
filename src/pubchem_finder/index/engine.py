"""Search the local index: FPSim2 for fingerprints, DuckDB for metadata,
RDKit to verify substructure candidates."""

import contextlib
import logging
import os
import sys
import threading
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Literal

import duckdb
import numpy as np
import psutil
from FPSim2 import FPSim2Engine
from rdkit import Chem

from .. import hints
from ..config import INDEX_MODE_ENV
from ..core.exceptions import IndexNotBuilt
from .layout import (
    ChunkStats,
    IndexPaths,
    Manifest,
    Merged,
    chunk_of,
    duckdb_path,
    unmerged,
)

logger = logging.getLogger(__name__)

Mode = Literal["auto", "memory", "disk"]
Metric = Literal["tanimoto", "dice", "cosine"]
# FPSim2 0.7.4 splits an in-memory search across threads, and its C++ search updates
# the fingerprint array's reference count from those threads unsafely: the count
# drifts, and when it drifts down the array is freed mid-search ("Chunk N worker
# died: ... GenericSearch(): incompatible function arguments" with a string where
# the array was, or a crash). FPSim2 runs the search in the calling thread when
# n_workers is 1, so that is the default. That is not enough on its own: the binding
# also drops its reference to the array without the GIL, so two searches of one
# in-memory array from different threads (parallel MCP tool calls) race the same way.
# LocalIndex therefore lets one in-memory search per array run at a time.
SEARCH_WORKERS = 1
# How long a chunk's stats file is read afresh on every search (see _chunk_stats):
# longer than any filesystem's modification-time step (FAT's is 2 s).
STATS_SETTLE_NS = 3_000_000_000
METADATA_COLUMNS = (
    "cid",
    "smiles",
    "inchikey",
    "iupac_name",
    "molecular_formula",
    "molecular_weight",
)


def available_memory() -> int | None:
    """Bytes of RAM available without swapping (Linux, macOS, Windows), or None
    if the platform does not say."""
    try:
        return int(psutil.virtual_memory().available)
    except (OSError, RuntimeError, AttributeError):
        return None


def choose_in_memory(needed_bytes: int, mode: Mode) -> bool:
    """In-memory search is much faster; use it when the fingerprints fit
    comfortably (with 25% headroom) in available RAM."""
    if mode != "auto":
        return mode == "memory"
    free = available_memory()
    return free is None or needed_bytes * 1.25 < free


def _open_engine(path: Path, in_memory: bool) -> Any:
    # FPSim2 print()s a warning when the file was built with another RDKit
    # version; stdout is the MCP transport, so send it to stderr instead.
    with contextlib.redirect_stdout(sys.stderr):
        return FPSim2Engine(str(path), in_memory_fps=in_memory)


@dataclass
class SubstructureMatches:
    """Verified substructure matches, in CID order, plus how far the search got."""

    merged: Merged
    """The merge whose pattern fingerprints gave the candidates."""
    rows: list[dict[str, Any]] = field(default_factory=list)
    n_candidates: int = 0
    n_checked: int = 0
    n_unverified: int = 0
    """Candidates with no stored structure to check: their chunk was rebuilt, or
    removed, after the merge. They are counted in ``n_checked``."""
    exhaustive: bool = True


class LocalIndex:
    """
    A finalized local index: similarity search, and substructure search when
    the index was built with pattern fingerprints.

    Each FPSim2 engine loads lazily on its first search (loading the full PubChem
    index takes tens of seconds and ~20 GB of RAM in memory mode), under a lock
    so concurrent callers load it once. In-memory searches run on one thread and one
    at a time per engine (see SEARCH_WORKERS); on-disk searches read fresh arrays
    and need neither.
    """

    def __init__(
        self, home: Path, mode: Mode | None = None, n_workers: int | None = None
    ):
        self.paths = IndexPaths(home)
        env_mode = os.environ.get(INDEX_MODE_ENV, "auto")
        self.mode: Mode = mode or (
            env_mode if env_mode in ("memory", "disk") else "auto"
        )  # type: ignore[assignment]
        self.n_workers = n_workers or SEARCH_WORKERS
        self._lock = threading.Lock()
        # One per engine, so a long substructure screen doesn't hold up similarity
        # searches; separate from _lock, so a search doesn't block loading.
        self._search_lock = threading.Lock()
        self._pattern_search_lock = threading.Lock()
        self._engine: Any = None
        self._in_memory: bool | None = None
        self._merged: Merged | None = None
        self._pattern_engine: Any = None
        self._pattern_in_memory: bool | None = None
        self._pattern_merged: Merged | None = None
        self._db: duckdb.DuckDBPyConnection | None = None
        self._db_files: list[Path] = []
        # Chunk stats already read, by chunk name, with what os.stat said then.
        self._stats_seen: dict[str, tuple[tuple[int, int, int], ChunkStats]] = {}
        # Its own lock: _lock is held while an engine loads, which can take long.
        self._stats_lock = threading.Lock()

    @property
    def merged(self) -> Merged:
        """The manifest on disk now, and when it was written."""
        merged = self.paths.read_manifest()
        if merged is None:
            raise IndexNotBuilt(
                f"No local index in {self.paths.home}. "
                + hints.build_index(self.paths.home)
            )
        return merged

    @property
    def manifest(self) -> Manifest:
        return self.merged.manifest

    @property
    def loaded(self) -> bool:
        return self._engine is not None

    @property
    def in_memory(self) -> bool | None:
        return self._in_memory

    def _ensure_loaded(self) -> tuple[Any, Merged]:
        """The similarity engine, and the merge it was loaded under."""
        with self._lock:
            merged = self.merged
            if self._engine is None or self._merged != merged:
                in_memory = choose_in_memory(merged.manifest.memory_bytes(), self.mode)
                logger.info(
                    "loading %s (%s mode)",
                    self.paths.fingerprints,
                    "memory" if in_memory else "disk",
                )
                self._engine = _open_engine(self.paths.fingerprints, in_memory)
                self._in_memory = in_memory
                self._merged = merged
            return self._engine, merged

    def search(
        self, query: Chem.Mol, threshold: float, top_k: int, metric: Metric = "tanimoto"
    ) -> list[tuple[int, float]]:
        """Top ``top_k`` (cid, similarity) pairs at or above ``threshold``, best first.

        Each CID appears once. Chunks can overlap (an ``--sdf`` file repeating
        compounds of a PubChem chunk), and the copies of a CID take FPSim2's slots,
        so the search asks for more until it holds ``top_k`` distinct CIDs and the
        tie at the cut is complete, or no more hits exist; without duplicates and
        without a tie at the cut that is a single FPSim2 call.
        """
        return self.search_merged(query, threshold, top_k, metric)[0]

    def search_merged(
        self, query: Chem.Mol, threshold: float, top_k: int, metric: Metric = "tanimoto"
    ) -> tuple[list[tuple[int, float]], Merged]:
        """:meth:`search`, and the merge that answered it.

        Pass that on to :meth:`metadata` and :meth:`unmerged`. A build can finish
        between the search and the lookup of its hits; judged against the manifest
        then on disk, a hit from the old fingerprints would be given the structure
        the new file has under the same ID.
        """
        engine, merged = self._ensure_loaded()
        # One past the cut: ask for top_k + 1 so the score after the cut is in
        # hand. If it ties with the score at the cut, compounds FPSim2 left out
        # could still win on CID; if overlapping chunks leave fewer than top_k
        # distinct CIDs, copies took slots. Either way, double the ask until the
        # lowest score in hand is below the score at the cut (the tie is complete,
        # so nothing omitted can reach it) or the hits run out (issue #159: the
        # cut used to follow FPSim2's own order).
        k = top_k + 1
        while True:
            hits = self._top_k(engine, query, threshold, k, metric)
            best: dict[int, float] = {}
            for hit in hits:
                cid, sim = int(hit["mol_id"]), round(float(hit["coeff"]), 4)
                best[cid] = max(sim, best.get(cid, sim))
            if len(hits) < k:
                break
            if len(best) > top_k:
                scores = sorted(best.values(), reverse=True)
                if scores[-1] < scores[top_k - 1]:
                    break
            k *= 2
        # FPSim2 returns equal scores in whatever order its threads finish. Break
        # ties by CID so results are stable and long-known compounds come first
        # (e.g. aspirin, CID 2244, before its 11C-labelled twin, CID 450661).
        ranked = sorted(best.items(), key=lambda pair: (-pair[1], pair[0]))[:top_k]
        return ranked, merged

    def _top_k(
        self, engine: Any, query: Chem.Mol, threshold: float, k: int, metric: Metric
    ) -> Any:
        if self._in_memory:
            with self._search_lock:
                return engine.top_k(
                    query,
                    k=k,
                    threshold=threshold,
                    metric=metric,
                    n_workers=self.n_workers,
                )
        # FPSim2 parallelizes on-disk search with forked processes, which can
        # deadlock when called from a threaded server (the MCP server runs
        # searches on worker threads), so on-disk search stays single-process.
        return engine.on_disk_top_k(
            query, k=k, threshold=threshold, metric=metric, n_workers=1
        )

    def _ensure_pattern_loaded(self) -> tuple[Any, Merged]:
        with self._lock:
            merged = self.merged
            sub = merged.manifest.substructure
            if sub is None or not self.paths.pattern_fingerprints.exists():
                raise IndexNotBuilt(
                    f"The local index in {self.paths.home} has no substructure search. "
                    f"Add it with {hints.add_substructure(self.paths.home)}."
                )
            if self._pattern_engine is None or self._pattern_merged != merged:
                in_memory = choose_in_memory(sub.memory_bytes(), self.mode)
                logger.info(
                    "loading %s (%s mode)",
                    self.paths.pattern_fingerprints,
                    "memory" if in_memory else "disk",
                )
                self._pattern_engine = _open_engine(
                    self.paths.pattern_fingerprints, in_memory
                )
                self._pattern_in_memory = in_memory
                self._pattern_merged = merged
            return self._pattern_engine, merged

    def screen(self, query: Chem.Mol) -> np.ndarray:
        """Sorted, distinct CIDs whose pattern fingerprint has every bit of the
        query's: a superset of the compounds containing ``query``. Distinct because
        overlapping chunks can hold a CID more than once."""
        return self._screen(query)[0]

    def _screen(self, query: Chem.Mol) -> tuple[np.ndarray, Merged]:
        engine, merged = self._ensure_pattern_loaded()
        # Fingerprint the query here rather than letting FPSim2 parse it, so
        # SMARTS query molecules work too.
        pattern_fp = Chem.PatternFingerprint(query)
        if self._pattern_in_memory:
            with self._pattern_search_lock:
                cids = engine.substructure(pattern_fp, n_workers=self.n_workers)
        else:  # single process on disk, see search()
            cids = engine.on_disk_substructure(pattern_fp, n_workers=1)
        return np.unique(np.asarray(cids, dtype=np.int64)), merged

    def substructure(
        self,
        query: Chem.Mol,
        max_hits: int,
        max_candidates: int,
        max_batch: int = 5000,
    ) -> SubstructureMatches:
        """
        Compounds containing ``query``, lowest CID first.

        Candidates from :meth:`screen` are verified with RDKit against their stored
        SMILES, in CID order, until ``max_hits`` matches are found or
        ``max_candidates`` have been checked. Metadata is fetched in batches that
        start small (most queries finish early) and double up to ``max_batch``.
        """
        candidates, merged = self._screen(query)
        result = SubstructureMatches(merged, n_candidates=int(candidates.size))
        limit = min(candidates.size, max_candidates)
        size = min(max_batch, max(256, 2 * max_hits))
        while result.n_checked < limit and len(result.rows) < max_hits:
            start = result.n_checked
            batch = [int(c) for c in candidates[start : min(start + size, limit)]]
            rows = self.metadata(batch, merged)
            for cid in batch:
                result.n_checked += 1
                row = rows.get(cid)
                if row is None or not row["smiles"]:
                    result.n_unverified += 1
                    continue
                mol = Chem.MolFromSmiles(row["smiles"])
                if mol is not None and mol.HasSubstructMatch(query):
                    result.rows.append(row)
                    if len(result.rows) >= max_hits:
                        break
            size = min(size * 2, max_batch)
        # A candidate that could not be checked may be a match: not the whole answer.
        result.exhaustive = (
            result.n_checked == result.n_candidates and not result.n_unverified
        )
        return result

    def unmerged(self, merged: Merged | None = None) -> dict[str, ChunkStats | None]:
        """Chunks rebuilt, or removed, since ``merged`` (default: the index on disk
        now; see :func:`~pubchem_finder.index.layout.unmerged`), each with its
        stats, or None when it is not complete on disk.

        The index still holds their earlier fingerprints, so :meth:`metadata` leaves
        their Parquet files out: a hit from one comes back without a structure
        rather than with the structure another build filed under the same ID.
        """
        merged = merged or self.merged
        seen = self._chunk_stats(merged.manifest.chunks)
        stats = {name: found for name, (found, _) in seen.items()}
        written = {name: at for name, (_, at) in seen.items()}
        stale = unmerged(merged, stats.values(), written)
        return {name: stats.get(name) for name in stale}

    def _chunk_stats(self, names: list[str]) -> dict[str, tuple[ChunkStats, int]]:
        """The stats of the chunks in ``names`` that are complete on disk, each with
        the time its file was written.

        One os.stat per chunk and search; a stats file is parsed again only when it
        was replaced (a build writes a new file, so its inode changes with it), or
        while it is so new that a replacement could look the same: where inode
        numbers are reused and modification times are coarse, a rebuild within the
        same tick leaves a file with the same number, time and size.
        """
        found: dict[str, tuple[ChunkStats, int]] = {}
        with self._stats_lock:
            for name in names:
                path = self.paths.chunk_stats(name)
                try:
                    stat = path.stat()
                    key = (stat.st_ino, stat.st_mtime_ns, stat.st_size)
                    seen = self._stats_seen.get(name)
                    settled = time.time_ns() - stat.st_mtime_ns > STATS_SETTLE_NS
                    if seen is None or seen[0] != key or not settled:
                        text = path.read_text(encoding="utf-8")
                        seen = (key, ChunkStats.model_validate_json(text))
                        self._stats_seen[name] = seen
                except (OSError, ValueError):  # no stats: not a complete chunk
                    self._stats_seen.pop(name, None)
                    continue
                found[name] = (seen[1], seen[0][1])
        return found

    def metadata(
        self, cids: list[int], merged: Merged | None = None
    ) -> dict[int, dict[str, Any]]:
        """Metadata rows for ``cids`` from the chunk Parquet files, leaving out
        chunks rebuilt since ``merged`` (see :meth:`unmerged`).

        Give the merge the IDs came from (:meth:`search_merged`); the default is
        the index on disk now. Every row returned was read from a chunk that was
        the merged build both before and after the read: a row from a chunk that
        a build replaced meanwhile is dropped, however many chunks are changing.
        """
        if not cids:
            return {}
        merged = merged or self.merged
        for attempt in range(3):
            stale = set(self.unmerged(merged))
            chunks = [n for n in merged.manifest.chunks if n not in stale]
            try:
                found = self._metadata(cids, chunks)
            except duckdb.Error:
                # A build replaced or removed a file while it was being read.
                if attempt == 2:
                    raise
                with self._lock:
                    self._db = None
                continue
            stale |= set(self.unmerged(merged))
            rows: dict[int, dict[str, Any]] = {}
            for chunk, row in found:
                if chunk in stale:
                    continue
                # One row per ID; one with a structure before one without.
                kept = rows.setdefault(int(row["cid"]), row)
                if not kept["smiles"] and row["smiles"]:
                    rows[int(row["cid"])] = row
            return rows
        return {}  # not reached: the third failure raises

    def _metadata(
        self, cids: list[int], chunks: list[str]
    ) -> list[tuple[str, dict[str, Any]]]:
        """Rows for ``cids``, each with the chunk whose file it was read from."""
        files = self.paths.metadata_files(chunks)
        if not files:
            return []
        with self._lock:
            if self._db is None or self._db_files != files:
                self._db = duckdb.connect(":memory:")
                self._db_files = files
                # A relation, not SQL text, so file paths need no quoting. By name:
                # chunks written by different versions may not share every column.
                self._db.read_parquet(
                    [duckdb_path(f) for f in files], union_by_name=True, filename=True
                ).create_view("meta")
            cursor = self._db.cursor()
        ids = ", ".join(str(int(c)) for c in cids)
        # Integers only; a literal IN list lets DuckDB push the filter into the scan.
        rows = cursor.execute(
            f"SELECT {', '.join(METADATA_COLUMNS)}, filename FROM meta "  # noqa: S608
            f"WHERE cid IN ({ids})"
        ).fetchall()
        return [
            (chunk_of(r[-1]), dict(zip(METADATA_COLUMNS, r[:-1], strict=True)))
            for r in rows
        ]
