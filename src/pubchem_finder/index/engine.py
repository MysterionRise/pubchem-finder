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

from ..config import INDEX_MODE_ENV
from ..core.exceptions import IndexNotBuilt
from .layout import IndexPaths, Manifest

logger = logging.getLogger(__name__)

Mode = Literal["auto", "memory", "disk"]
Metric = Literal["tanimoto", "dice", "cosine"]
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
    so concurrent callers load it once.
    """

    def __init__(
        self, home: Path, mode: Mode | None = None, n_workers: int | None = None
    ):
        self.paths = IndexPaths(home)
        env_mode = os.environ.get(INDEX_MODE_ENV, "auto")
        self.mode: Mode = mode or (
            env_mode if env_mode in ("memory", "disk") else "auto"
        )  # type: ignore[assignment]
        self.n_workers = n_workers or max(1, min(4, (os.cpu_count() or 1)))
        self._lock = threading.Lock()
        self._engine: Any = None
        self._in_memory: bool | None = None
        self._manifest: Manifest | None = None
        self._pattern_engine: Any = None
        self._pattern_in_memory: bool | None = None
        self._pattern_manifest: Manifest | None = None
        self._db: duckdb.DuckDBPyConnection | None = None
        self._db_files: list[Path] = []

    @property
    def manifest(self) -> Manifest:
        manifest = self.paths.load_manifest()
        if manifest is None:
            raise IndexNotBuilt(
                f"No local index in {self.paths.home}. Build one with "
                "`pubchem-finder index build --limit 1` (a ~500k-compound sample, minutes) "
                "or `pubchem-finder index build --all` (all of PubChem, an overnight job)."
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
        """Top ``top_k`` (cid, similarity) pairs at or above ``threshold``, best first."""
        engine = self._ensure_loaded()
        if self._in_memory:
            hits = engine.top_k(
                query,
                k=top_k,
                threshold=threshold,
                metric=metric,
                n_workers=self.n_workers,
            )
        else:
            # FPSim2 parallelizes on-disk search with forked processes, which can
            # deadlock when called from a threaded server (the MCP server runs
            # searches on worker threads), so on-disk search stays single-process.
            hits = engine.on_disk_top_k(
                query, k=top_k, threshold=threshold, metric=metric, n_workers=1
            )
        return [(int(h["mol_id"]), round(float(h["coeff"]), 4)) for h in hits]

    def _ensure_pattern_loaded(self) -> Any:
        with self._lock:
            manifest = self.manifest
            sub = manifest.substructure
            if sub is None or not self.paths.pattern_fingerprints.exists():
                raise IndexNotBuilt(
                    f"The local index in {self.paths.home} has no substructure search. "
                    "Add it with `pubchem-finder index build --with-substructure` "
                    "(uses the compounds already indexed; no download)."
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
        """Sorted CIDs whose pattern fingerprint has every bit of the query's: a
        superset of the compounds containing ``query``."""
        engine = self._ensure_pattern_loaded()
        # Fingerprint the query here rather than letting FPSim2 parse it, so
        # SMARTS query molecules work too.
        pattern_fp = Chem.PatternFingerprint(query)
        if self._pattern_in_memory:
            cids = engine.substructure(pattern_fp, n_workers=self.n_workers)
        else:  # single process on disk, see search()
            cids = engine.on_disk_substructure(pattern_fp, n_workers=1)
        return np.sort(np.asarray(cids, dtype=np.int64))

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

    def metadata(self, cids: list[int]) -> dict[int, dict[str, Any]]:
        """Metadata rows for ``cids`` from the chunk Parquet files."""
        if not cids:
            return {}
        manifest = self.manifest
        files = self.paths.metadata_files(manifest.chunks)
        if not files:
            return {}
        with self._lock:
            if self._db is None or self._db_files != files:
                self._db = duckdb.connect(":memory:")
                self._db_files = files
                file_list = ", ".join(
                    "'" + str(f).replace("'", "''") + "'" for f in files
                )
                self._db.execute(
                    f"CREATE VIEW meta AS SELECT * FROM read_parquet([{file_list}])"
                )
            cursor = self._db.cursor()
        ids = ", ".join(str(int(c)) for c in cids)
        rows = cursor.execute(
            f"SELECT {', '.join(METADATA_COLUMNS)} FROM meta WHERE cid IN ({ids})"
        ).fetchall()
        return {int(r[0]): dict(zip(METADATA_COLUMNS, r, strict=True)) for r in rows}
