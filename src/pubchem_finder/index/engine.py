"""Search the local index: FPSim2 for fingerprints, DuckDB for metadata,
RDKit to verify substructure candidates."""

import contextlib
import logging
import os
import sys
import threading
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
from .layout import ChunkStats, IndexPaths, Manifest, unmerged

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

    rows: list[dict[str, Any]] = field(default_factory=list)
    n_candidates: int = 0
    n_checked: int = 0
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
        self._manifest: Manifest | None = None
        self._pattern_engine: Any = None
        self._pattern_in_memory: bool | None = None
        self._pattern_manifest: Manifest | None = None
        self._db: duckdb.DuckDBPyConnection | None = None
        self._db_files: list[Path] = []
        # Chunk stats already read, by chunk name, with what os.stat said then.
        self._stats_seen: dict[str, tuple[tuple[int, int, int], ChunkStats]] = {}
        # Its own lock: _lock is held while an engine loads, which can take long.
        self._stats_lock = threading.Lock()

    @property
    def manifest(self) -> Manifest:
        manifest = self.paths.load_manifest()
        if manifest is None:
            raise IndexNotBuilt(
                f"No local index in {self.paths.home}. "
                + hints.build_index(self.paths.home)
            )
        return manifest

    @property
    def loaded(self) -> bool:
        return self._engine is not None

    @property
    def in_memory(self) -> bool | None:
        return self._in_memory

    def _ensure_loaded(self) -> Any:
        with self._lock:
            manifest = self.manifest
            if self._engine is None or self._manifest != manifest:
                in_memory = choose_in_memory(manifest.memory_bytes(), self.mode)
                logger.info(
                    "loading %s (%s mode)",
                    self.paths.fingerprints,
                    "memory" if in_memory else "disk",
                )
                self._engine = _open_engine(self.paths.fingerprints, in_memory)
                self._in_memory = in_memory
                self._manifest = manifest
            return self._engine

    def search(
        self, query: Chem.Mol, threshold: float, top_k: int, metric: Metric = "tanimoto"
    ) -> list[tuple[int, float]]:
        """Top ``top_k`` (cid, similarity) pairs at or above ``threshold``, best first.

        Each CID appears once. Chunks can overlap (an ``--sdf`` file repeating
        compounds of a PubChem chunk), and the copies of a CID take FPSim2's slots,
        so the search asks for more until it has ``top_k`` distinct CIDs or no
        more hits exist; without duplicates that is a single FPSim2 call.
        """
        engine = self._ensure_loaded()
        k = top_k
        while True:
            hits = self._top_k(engine, query, threshold, k, metric)
            best: dict[int, float] = {}
            for hit in hits:
                cid, sim = int(hit["mol_id"]), round(float(hit["coeff"]), 4)
                best[cid] = max(sim, best.get(cid, sim))
            if len(best) >= top_k or len(hits) < k:
                break
            k *= 2
        # FPSim2 returns equal scores in whatever order its threads finish. Break
        # ties by CID so results are stable and long-known compounds come first
        # (e.g. aspirin, CID 2244, before its 11C-labelled twin, CID 450661).
        return sorted(best.items(), key=lambda pair: (-pair[1], pair[0]))[:top_k]

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

    def _ensure_pattern_loaded(self) -> Any:
        with self._lock:
            manifest = self.manifest
            sub = manifest.substructure
            if sub is None or not self.paths.pattern_fingerprints.exists():
                raise IndexNotBuilt(
                    f"The local index in {self.paths.home} has no substructure search. "
                    f"Add it with {hints.add_substructure(self.paths.home)}."
                )
            if self._pattern_engine is None or self._pattern_manifest != manifest:
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
                self._pattern_manifest = manifest
            return self._pattern_engine

    def screen(self, query: Chem.Mol) -> np.ndarray:
        """Sorted, distinct CIDs whose pattern fingerprint has every bit of the
        query's: a superset of the compounds containing ``query``. Distinct because
        overlapping chunks can hold a CID more than once."""
        engine = self._ensure_pattern_loaded()
        # Fingerprint the query here rather than letting FPSim2 parse it, so
        # SMARTS query molecules work too.
        pattern_fp = Chem.PatternFingerprint(query)
        if self._pattern_in_memory:
            with self._pattern_search_lock:
                cids = engine.substructure(pattern_fp, n_workers=self.n_workers)
        else:  # single process on disk, see search()
            cids = engine.on_disk_substructure(pattern_fp, n_workers=1)
        return np.unique(np.asarray(cids, dtype=np.int64))

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
        candidates = self.screen(query)
        result = SubstructureMatches(n_candidates=int(candidates.size))
        limit = min(candidates.size, max_candidates)
        size = min(max_batch, max(256, 2 * max_hits))
        while result.n_checked < limit and len(result.rows) < max_hits:
            start = result.n_checked
            batch = [int(c) for c in candidates[start : min(start + size, limit)]]
            rows = self.metadata(batch)
            for cid in batch:
                result.n_checked += 1
                row = rows.get(cid)
                if row is None or not row["smiles"]:
                    continue
                mol = Chem.MolFromSmiles(row["smiles"])
                if mol is not None and mol.HasSubstructMatch(query):
                    result.rows.append(row)
                    if len(result.rows) >= max_hits:
                        break
            size = min(size * 2, max_batch)
        result.exhaustive = result.n_checked == result.n_candidates
        return result

    def unmerged(self, manifest: Manifest | None = None) -> dict[str, str | None]:
        """Chunks rebuilt, or removed, since the index was merged (see
        :func:`~pubchem_finder.index.layout.unmerged`), each with the file it was
        last built from when that is known.

        The index still holds their earlier fingerprints, so :meth:`metadata` leaves
        their Parquet files out: a hit from one comes back without a structure
        rather than with the structure another build filed under the same ID.
        """
        manifest = manifest or self.manifest
        stats = {c.name: c for c in self._chunk_stats(manifest.chunks)}
        files: dict[str, str | None] = {}
        for name in unmerged(manifest, stats.values()):
            source = stats[name].source if name in stats else None
            files[name] = source.path if source is not None else None
        return files

    def _chunk_stats(self, names: list[str]) -> list[ChunkStats]:
        """The stats of the chunks in ``names`` that are complete on disk.

        One os.stat per chunk and search; a stats file is parsed again only when it
        was replaced (a build writes a new file, so its inode changes with it).
        """
        found: list[ChunkStats] = []
        with self._stats_lock:
            for name in names:
                path = self.paths.chunk_stats(name)
                try:
                    stat = path.stat()
                    key = (stat.st_ino, stat.st_mtime_ns, stat.st_size)
                    seen = self._stats_seen.get(name)
                    if seen is None or seen[0] != key:
                        text = path.read_text(encoding="utf-8")
                        seen = (key, ChunkStats.model_validate_json(text))
                        self._stats_seen[name] = seen
                except (OSError, ValueError):  # no stats: not a complete chunk
                    self._stats_seen.pop(name, None)
                    continue
                found.append(seen[1])
        return found

    def metadata(self, cids: list[int]) -> dict[int, dict[str, Any]]:
        """Metadata rows for ``cids`` from the chunk Parquet files, leaving out
        chunks rebuilt since the index was merged (see :meth:`unmerged`)."""
        if not cids:
            return {}
        rows: dict[int, dict[str, Any]] = {}
        for _ in range(3):
            manifest = self.manifest
            stale = set(self.unmerged(manifest))
            rows = self._metadata(cids, [n for n in manifest.chunks if n not in stale])
            # A build may have replaced a chunk while its file was being read.
            if set(self.unmerged(manifest)) == stale:
                break
        return rows

    def _metadata(
        self, cids: list[int], chunks: list[str]
    ) -> dict[int, dict[str, Any]]:
        files = self.paths.metadata_files(chunks)
        if not files:
            return {}
        with self._lock:
            if self._db is None or self._db_files != files:
                self._db = duckdb.connect(":memory:")
                self._db_files = files
                # A relation, not SQL text, so file paths need no quoting. By name:
                # chunks written by different versions may not share every column.
                self._db.read_parquet(
                    [str(f) for f in files], union_by_name=True
                ).create_view("meta")
            cursor = self._db.cursor()
        ids = ", ".join(str(int(c)) for c in cids)
        # Integers only; a literal IN list lets DuckDB push the filter into the scan.
        rows = cursor.execute(
            f"SELECT {', '.join(METADATA_COLUMNS)} FROM meta WHERE cid IN ({ids})"  # noqa: S608
        ).fetchall()
        return {int(r[0]): dict(zip(METADATA_COLUMNS, r, strict=True)) for r in rows}
