"""On-disk cache for PubChem web-service responses.

PubChem data changes slowly, and an assistant tends to ask the same questions
several times in one conversation. Caching successful responses for a week
makes repeats free and keeps the client well inside PubChem's request limits.
Only successful (HTTP 200) bodies are stored; errors are never cached.

The cache is a single SQLite file in WAL mode, so the CLI and a running MCP
server can share it safely. It stays bounded: expired entries are dropped when
it is opened and now and then as it is written, and once the stored bodies pass
a size limit (256 MB by default, ``PUBCHEM_FINDER_CACHE_MAX_MB``) the oldest
entries go first.
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


@dataclass(frozen=True)
class CacheStats:
    path: Path
    entries: int
    size_bytes: int
    max_bytes: int


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
        max_bytes: int = DEFAULT_MAX_BYTES,
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
            self._maintain(con)

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
        with self._connect() as con:
            con.execute(
                "INSERT OR REPLACE INTO responses (key, body, created) VALUES (?, ?, ?)",
                (key, zlib.compress(body), self._clock()),
            )
            self._writes_since_maintenance += 1
            if self._writes_since_maintenance >= _MAINTAIN_EVERY:
                self._maintain(con)

    def _maintain(self, con: sqlite3.Connection) -> None:
        """Drop expired entries, then the oldest ones until under ``max_bytes``."""
        self._writes_since_maintenance = 0
        con.execute(
            "DELETE FROM responses WHERE created < ?",
            (self._clock() - self.ttl_seconds,),
        )
        rows = con.execute(
            "SELECT key, LENGTH(body) FROM responses ORDER BY created"
        ).fetchall()
        total = sum(size for _, size in rows)
        drop: list[tuple[str]] = []
        for key, size in rows:
            if total <= self.max_bytes:
                break
            drop.append((key,))
            total -= size
        con.executemany("DELETE FROM responses WHERE key = ?", drop)

    def clear(self) -> int:
        """Delete every entry; returns how many were removed."""
        with self._connect() as con:
            removed = con.execute("DELETE FROM responses").rowcount
        with closing(sqlite3.connect(self.path, timeout=10)) as con:
            con.execute("VACUUM")
        return int(removed)

    def stats(self) -> CacheStats:
        with self._connect() as con:
            (entries,) = con.execute("SELECT COUNT(*) FROM responses").fetchone()
        size = sum(
            p.stat().st_size
            for p in self.path.parent.glob(self.path.name + "*")
            if p.is_file()
        )
        return CacheStats(
            path=self.path,
            entries=int(entries),
            size_bytes=size,
            max_bytes=self.max_bytes,
        )
