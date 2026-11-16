"""On-disk cache for PubChem web-service responses.

PubChem data changes slowly, and an assistant tends to ask the same questions
several times in one conversation. Caching successful responses for a week
makes repeats free and keeps the client well inside PubChem's request limits.
Only successful (HTTP 200) bodies are stored; errors are never cached.

The cache is a single SQLite file in WAL mode, so the CLI and a running MCP
server can share it safely. It stays bounded: expired entries are dropped when
it is opened and now and then as it is written, and once the file passes a size
limit (256 MB by default, ``PUBCHEM_FINDER_CACHE_MAX_MB``) the oldest entries go
first. The limit is measured against the space the file uses, read straight from
SQLite's page counters, so the common case costs nothing; the row-by-row pass
runs only when the file is actually over the limit, and the freed space is
reclaimed with ``VACUUM`` so the file shrinks instead of staying at its largest.
"""

import hashlib
import json
import sqlite3
import threading
import time
import zlib
from collections.abc import Callable, Iterator, Mapping
from contextlib import closing, contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Any

DEFAULT_TTL_SECONDS = 7 * 24 * 3600
DEFAULT_MAX_BYTES = 256 * 1_000_000
# Prune expired entries and trim to size once every this many writes, so a
# long-lived cache stays bounded without paying for a sweep on each one.
_MAINTAIN_EVERY = 100
# Only VACUUM (which rewrites the whole file) when a maintenance pass freed at
# least this much, so routine pruning stays cheap but a big eviction shrinks
# the file instead of leaving it at its largest size.
_VACUUM_MIN_FREE_BYTES = 16 * 1_000_000


@dataclass(frozen=True)
class CacheStats:
    path: Path
    entries: int
    size_bytes: int
    max_bytes: float


class ResponseCache:
    """
    Args:
        path: SQLite file; parent directories are created.
        ttl_seconds: Entries older than this are ignored and overwritten.
        max_bytes: Soft cap on the stored bodies; the oldest entries are dropped
            when it is exceeded.
        clock: Time source (tests).
    """

    def __init__(
        self,
        path: Path,
        ttl_seconds: float = DEFAULT_TTL_SECONDS,
        max_bytes: float = DEFAULT_MAX_BYTES,
        clock: Callable[[], float] = time.time,
    ) -> None:
        self.path = Path(path)
        self.ttl_seconds = ttl_seconds
        self.max_bytes = max_bytes
        self._clock = clock
        self._lock = threading.Lock()
        self._writes_since_maintenance = 0
        self.path.parent.mkdir(parents=True, exist_ok=True)
        with self._connect() as con:
            con.execute("PRAGMA journal_mode=WAL")
            con.execute(
                "CREATE TABLE IF NOT EXISTS responses ("
                " key TEXT PRIMARY KEY, body BLOB NOT NULL, created REAL NOT NULL)"
            )
            # Index the column the expiry delete and the age-ordered trim both
            # use, so neither has to scan and sort the whole table.
            con.execute(
                "CREATE INDEX IF NOT EXISTS responses_created ON responses (created)"
            )
            reclaim = self._maintain(con)
        if reclaim:
            self._vacuum()

    @contextmanager
    def _connect(self) -> Iterator[sqlite3.Connection]:
        with self._lock, closing(sqlite3.connect(self.path, timeout=10)) as con, con:
            yield con

    @staticmethod
    def key(
        method: str,
        url: str,
        params: Mapping[str, Any] | None = None,
        data: Mapping[str, Any] | None = None,
    ) -> str:
        """Stable key for a request; parameter order does not matter."""
        canonical = json.dumps(
            [
                method.upper(),
                url,
                sorted((params or {}).items()),
                sorted((data or {}).items()),
            ],
            default=str,
        )
        return hashlib.sha256(canonical.encode()).hexdigest()

    def get(self, key: str) -> bytes | None:
        with self._connect() as con:
            row = con.execute(
                "SELECT body, created FROM responses WHERE key = ?", (key,)
            ).fetchone()
        if row is None or self._clock() - row[1] > self.ttl_seconds:
            return None
        return zlib.decompress(row[0])

    def put(self, key: str, body: bytes) -> None:
        reclaim = False
        with self._connect() as con:
            con.execute(
                "INSERT OR REPLACE INTO responses (key, body, created) VALUES (?, ?, ?)",
                (key, zlib.compress(body), self._clock()),
            )
            self._writes_since_maintenance += 1
            if self._writes_since_maintenance >= _MAINTAIN_EVERY:
                reclaim = self._maintain(con)
        if reclaim:
            self._vacuum()

    @staticmethod
    def _used_bytes(con: sqlite3.Connection) -> int:
        """Space the file actually holds, from SQLite's own page counters.

        ``(page_count - freelist_count) * page_size`` is the live data; it costs
        a header read instead of scanning the table, and freed-but-not-yet-
        reclaimed pages (the freelist) do not count against the limit.
        """
        page_size = int(con.execute("PRAGMA page_size").fetchone()[0])
        page_count = int(con.execute("PRAGMA page_count").fetchone()[0])
        freelist = int(con.execute("PRAGMA freelist_count").fetchone()[0])
        return (page_count - freelist) * page_size

    def _maintain(self, con: sqlite3.Connection) -> bool:
        """Drop expired entries, then the oldest ones until under ``max_bytes``.

        Returns whether enough was freed to be worth a ``VACUUM`` afterwards
        (which the caller runs, as it cannot run inside this transaction).
        """
        self._writes_since_maintenance = 0
        freed = con.execute(
            "DELETE FROM responses WHERE created < ?",
            (self._clock() - self.ttl_seconds,),
        ).rowcount
        # The cheap check: scan and sort the table only when the file is over
        # the limit, never on the common in-bounds open or write.
        if self._used_bytes(con) > self.max_bytes:
            freed += self._trim_to_size(con)
        if freed == 0:
            return False
        page_size = int(con.execute("PRAGMA page_size").fetchone()[0])
        freelist = int(con.execute("PRAGMA freelist_count").fetchone()[0])
        return freelist * page_size >= _VACUUM_MIN_FREE_BYTES

    def _trim_to_size(self, con: sqlite3.Connection) -> int:
        """Drop the oldest entries until the file is back under the limit."""
        rows = con.execute(
            "SELECT key, LENGTH(body) FROM responses ORDER BY created"
        ).fetchall()
        stored = sum(size for _, size in rows)
        if stored <= 0:
            return 0
        # The limit is on the file, which is larger than the stored bodies
        # (page and index overhead). Shed the same proportion of body bytes so
        # the file lands near the limit without measuring it row by row.
        target = stored * (self.max_bytes / self._used_bytes(con))
        drop: list[tuple[str]] = []
        for key, size in rows:
            if stored <= target:
                break
            drop.append((key,))
            stored -= size
        con.executemany("DELETE FROM responses WHERE key = ?", drop)
        return len(drop)

    def _vacuum(self) -> None:
        """Reclaim freed pages so the file shrinks. Runs outside a transaction."""
        with self._lock, closing(sqlite3.connect(self.path, timeout=10)) as con:
            con.execute("VACUUM")

    def clear(self) -> int:
        """Delete every entry; returns how many were removed."""
        with self._connect() as con:
            removed = con.execute("DELETE FROM responses").rowcount
        self._vacuum()
        return int(removed)

    def stats(self) -> CacheStats:
        with self._connect() as con:
            (entries,) = con.execute("SELECT COUNT(*) FROM responses").fetchone()
            size = self._used_bytes(con)
        return CacheStats(
            path=self.path,
            entries=int(entries),
            size_bytes=size,
            max_bytes=self.max_bytes,
        )
