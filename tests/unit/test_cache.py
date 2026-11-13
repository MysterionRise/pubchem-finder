"""Tests for the on-disk PubChem response cache."""

import httpx
import pytest
import respx

from pubchem_finder.data.cache import ResponseCache
from pubchem_finder.data.pug_rest import CompoundNotFound, PubChemClient
from tests import pubchem_fixtures as fx


class Clock:
    def __init__(self):
        self.now = 1_000_000.0

    def __call__(self):
        return self.now


def test_put_get_and_ttl(tmp_path):
    clock = Clock()
    cache = ResponseCache(tmp_path / "c" / "cache.sqlite", ttl_seconds=60, clock=clock)
    cache.put("k", b'{"a": 1}')
    assert cache.get("k") == b'{"a": 1}'
    assert cache.get("missing") is None
    clock.now += 61
    assert cache.get("k") is None  # expired
    cache.put("k", b"new")
    assert cache.get("k") == b"new"


def test_key_is_stable_and_order_independent():
    a = ResponseCache.key("get", "https://x", {"b": 2, "a": 1}, None)
    b = ResponseCache.key("GET", "https://x", {"a": 1, "b": 2}, {})
    assert a == b
    assert a != ResponseCache.key("POST", "https://x", {"a": 1, "b": 2}, None)
    assert a != ResponseCache.key("GET", "https://x", {"a": 1, "b": 3}, None)


def test_stats_and_clear(tmp_path):
    cache = ResponseCache(tmp_path / "cache.sqlite")
    cache.put("a", b"x" * 1000)
    cache.put("b", b"y")
    stats = cache.stats()
    assert stats.entries == 2
    assert stats.size_bytes > 0
    assert cache.clear() == 2
    assert cache.stats().entries == 0


def test_shared_between_instances(tmp_path):
    ResponseCache(tmp_path / "cache.sqlite").put("k", b"v")
    assert ResponseCache(tmp_path / "cache.sqlite").get("k") == b"v"


def test_expired_entries_are_removed_when_opened(tmp_path):
    clock = Clock()
    path = tmp_path / "cache.sqlite"
    ResponseCache(path, ttl_seconds=60, clock=clock).put("k", b"v")
    clock.now += 61
    reopened = ResponseCache(path, ttl_seconds=60, clock=clock)
    assert reopened.stats().entries == 0  # pruned on open, not just ignored


def test_size_limit_drops_the_oldest_first(tmp_path):
    import os

    clock = Clock()
    path = tmp_path / "cache.sqlite"
    bodies = {f"k{i}": os.urandom(4_000) for i in range(6)}  # incompressible
    cache = ResponseCache(path, max_bytes=10_000, clock=clock)
    for key, body in bodies.items():
        clock.now += 1
        cache.put(key, body)
    # Reopening trims to the cap, keeping the newest entries.
    trimmed = ResponseCache(path, max_bytes=10_000, clock=clock)
    stats = trimmed.stats()
    assert stats.max_bytes == 10_000
    assert trimmed.get("k0") is None  # oldest evicted
    assert trimmed.get("k5") == bodies["k5"]  # newest kept
    assert stats.entries <= 3  # only a few 4 KB bodies fit under 10 KB


def test_maintenance_runs_while_the_cache_is_used(tmp_path):
    from pubchem_finder.data import cache as cache_mod

    clock = Clock()
    cache = ResponseCache(tmp_path / "cache.sqlite", ttl_seconds=60, clock=clock)
    cache.put("old", b"v")
    clock.now += 61  # "old" is now expired
    for i in range(cache_mod._MAINTAIN_EVERY):
        cache.put(f"k{i}", b"v")
    # The periodic sweep deleted the expired entry, not only hid it.
    assert cache.get("old") is None
    assert cache.stats().entries == cache_mod._MAINTAIN_EVERY


def test_cache_max_bytes_from_the_environment(monkeypatch):
    from pubchem_finder.config import cache_max_bytes
    from pubchem_finder.data.cache import DEFAULT_MAX_BYTES

    monkeypatch.delenv("PUBCHEM_FINDER_CACHE_MAX_MB", raising=False)
    assert cache_max_bytes() == DEFAULT_MAX_BYTES
    monkeypatch.setenv("PUBCHEM_FINDER_CACHE_MAX_MB", "10")
    assert cache_max_bytes() == 10_000_000
    monkeypatch.setenv("PUBCHEM_FINDER_CACHE_MAX_MB", "nonsense")
    assert cache_max_bytes() == DEFAULT_MAX_BYTES


@pytest.fixture
def cached_client(tmp_path, monkeypatch):
    monkeypatch.setattr("pubchem_finder.data.pug_rest.time.sleep", lambda s: None)
    with PubChemClient(
        min_interval=0, cache=ResponseCache(tmp_path / "cache.sqlite")
    ) as c:
        yield c


def test_client_serves_repeats_from_cache(cached_client):
    with respx.mock() as router:
        route = router.get(fx.REST + "/compound/cid/2244/synonyms/JSON").respond(
            200, json=fx.SYNONYMS
        )
        first = cached_client.compound_json(2244, "synonyms")
        second = cached_client.compound_json(2244, "synonyms")
    assert first == second == fx.SYNONYMS
    assert route.call_count == 1


def test_post_body_is_part_of_the_key(cached_client):
    url = fx.REST + "/compound/fastidentity/smiles/property/X/JSON"
    with respx.mock() as router:
        route = router.post(url__regex=r".*fastidentity.*").respond(
            200, json=fx.table(fx.ASPIRIN)
        )
        cached_client.request("POST", url, data={"smiles": "CCO"})
        cached_client.request("POST", url, data={"smiles": "CCO"})
        cached_client.request("POST", url, data={"smiles": "CCN"})
    assert route.call_count == 2


def test_errors_are_not_cached(cached_client):
    with respx.mock() as router:
        route = router.get(fx.REST + "/compound/cid/2244/synonyms/JSON")
        route.side_effect = [
            httpx.Response(404, json=fx.NOT_FOUND),
            httpx.Response(200, json=fx.SYNONYMS),
        ]
        with pytest.raises(CompoundNotFound):
            cached_client.compound_json(2244, "synonyms")
        assert cached_client.compound_json(2244, "synonyms") == fx.SYNONYMS
    assert route.call_count == 2


# --- A damaged or unwritable cache falls back to the network ------------------

SYN_URL = fx.REST + "/compound/cid/2244/synonyms/JSON"


def test_cache_failures_fall_back_to_network(cached_client, monkeypatch):
    import sqlite3

    def broken(*args):
        raise sqlite3.OperationalError("database is locked")

    monkeypatch.setattr(cached_client.cache, "get", broken)
    monkeypatch.setattr(cached_client.cache, "put", broken)
    with respx.mock() as router:
        route = router.get(SYN_URL).respond(200, json=fx.SYNONYMS)
        assert cached_client.compound_json(2244, "synonyms") == fx.SYNONYMS
        assert cached_client.compound_json(2244, "synonyms") == fx.SYNONYMS
    assert route.call_count == 2


def test_corrupt_cache_entry_is_refetched(cached_client):
    import sqlite3
    import zlib

    with respx.mock() as router:
        route = router.get(SYN_URL).respond(200, json=fx.SYNONYMS)
        cached_client.compound_json(2244, "synonyms")
        with sqlite3.connect(cached_client.cache.path) as con:
            con.execute("UPDATE responses SET body = ?", (zlib.compress(b"\xffjunk"),))
        assert cached_client.compound_json(2244, "synonyms") == fx.SYNONYMS
        assert cached_client.compound_json(2244, "synonyms") == fx.SYNONYMS
    assert route.call_count == 2  # refetched once, then served from the cache again


def test_damaged_cache_file_falls_back_to_network(cached_client):
    cached_client.cache.path.write_bytes(b"this is not a database" * 100)
    with respx.mock() as router:
        router.get(SYN_URL).respond(200, json=fx.SYNONYMS)
        assert cached_client.compound_json(2244, "synonyms") == fx.SYNONYMS


def test_only_json_objects_are_cached(cached_client):
    with respx.mock() as router:
        route = router.get(SYN_URL).respond(200, json=[])
        assert cached_client.request("GET", SYN_URL) == []
        assert cached_client.request("GET", SYN_URL) == []
    assert route.call_count == 2
