"""Search the local index: FPSim2 for fingerprints, DuckDB for metadata."""

import logging
import os
import threading
from pathlib import Path
from typing import Any, Literal

import duckdb
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
    """Best-effort bytes of available RAM, or None if unknown."""
    try:
        with open("/proc/meminfo") as fh:
            for line in fh:
                if line.startswith("MemAvailable:"):
                    return int(line.split()[1]) * 1024
    except OSError:
        pass
    try:
        return os.sysconf("SC_AVPHYS_PAGES") * os.sysconf("SC_PAGE_SIZE")
    except (AttributeError, OSError, ValueError):
        return None


def choose_in_memory(manifest: Manifest, mode: Mode) -> bool:
    """In-memory search is much faster; use it when the fingerprints fit
    comfortably (with 25% headroom) in available RAM."""
    if mode != "auto":
        return mode == "memory"
    free = available_memory()
    return free is None or manifest.memory_bytes() * 1.25 < free


class SimilarityIndex:
    """
    A finalized local index.

    The FPSim2 engine loads lazily on the first search (loading the full PubChem
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
        self._db: duckdb.DuckDBPyConnection | None = None

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
                in_memory = choose_in_memory(manifest, self.mode)
                logger.info(
                    "loading %s (%s mode)",
                    self.paths.fingerprints,
                    "memory" if in_memory else "disk",
                )
                self._engine = FPSim2Engine(
                    str(self.paths.fingerprints), in_memory_fps=in_memory
                )
                self._in_memory = in_memory
                self._manifest = manifest
                self._db = None
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

    def metadata(self, cids: list[int]) -> dict[int, dict[str, Any]]:
        """Metadata rows for ``cids`` from the chunk Parquet files."""
        if not cids:
            return {}
        manifest = self.manifest
        files = self.paths.metadata_files(manifest.chunks)
        if not files:
            return {}
        with self._lock:
            if self._db is None:
                self._db = duckdb.connect(":memory:")
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
