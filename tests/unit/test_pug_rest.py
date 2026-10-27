"""Tests for the PUG REST client (all HTTP mocked with respx)."""

import os
import time

import httpx
import pytest
import respx

from pubchem_finder.data.pug_rest import (
    PROPERTIES,
    BadQuery,
    CompoundNotFound,
    PubChemClient,
    PubChemError,
    PubChemUnavailable,
)
from tests.conftest import CAFFEINE_ROW, fault, property_table


def test_lookup_by_name(pubchem_api, fast_client):
    records = fast_client.lookup("caffeine")
    assert len(records) == 1
    rec = records[0]
    assert rec.cid == 2519
    assert rec.title == "Caffeine"
    assert rec.molecular_weight == pytest.approx(194.19)  # coerced from string
    assert rec.smiles == "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
    assert rec.inchikey == "RYYVLZVUVIJVGH-UHFFFAOYSA-N"
    assert rec.pubchem_url == "https://pubchem.ncbi.nlm.nih.gov/compound/2519"
    request = pubchem_api.calls.last.request
    assert request.method == "GET"
    assert "/compound/name/caffeine/property/" in str(request.url)
    assert ",".join(PROPERTIES) in str(request.url)
    assert request.headers["user-agent"].startswith("pubchem-finder/")


def test_lookup_by_cid_and_limit(pubchem_api, fast_client):
    records = fast_client.lookup("2244,2519", by="cid", limit=1)
    assert [r.cid for r in records] == [2244]
    records = fast_client.lookup("2244,2519", by="cid", limit=5)
    assert [r.cid for r in records] == [2244, 2519]


def test_lookup_by_inchikey(pubchem_api, fast_client):
    records = fast_client.lookup("BSYNRYMUTXBXSQ-UHFFFAOYSA-N", by="inchikey")
    assert records[0].cid == 2244


def test_lookup_by_smiles_uses_post_form(pubchem_api, fast_client):
    records = fast_client.lookup("CC(=O)Oc1ccccc1C(=O)O", by="smiles")
    assert records[0].cid == 2244
    request = pubchem_api.calls.last.request
    assert request.method == "POST"
    assert request.content == b"smiles=CC%28%3DO%29Oc1ccccc1C%28%3DO%29O"


def test_name_with_spaces_is_url_quoted(fast_client):
    with respx.mock(base_url="https://pubchem.ncbi.nlm.nih.gov/rest/pug") as router:
        route = router.get(
            path__regex=r"/compound/name/acetylsalicylic acid/property/.*"
        )
        route.respond(200, json=property_table(CAFFEINE_ROW))
        fast_client.lookup("acetylsalicylic acid")
        assert route.called
        assert "/name/acetylsalicylic%20acid/" in str(route.calls.last.request.url)


def test_not_found_raises(pubchem_api, fast_client):
    with pytest.raises(CompoundNotFound, match="No CID found"):
        fast_client.lookup("notacompound")


def test_bad_request_raises(pubchem_api, fast_client):
    with pytest.raises(BadQuery, match="Invalid CID"):
        fast_client.lookup("abc", by="cid")


@pytest.mark.parametrize(
    ("query", "by", "limit"), [("", "name", 5), ("x", "bogus", 5), ("x", "name", 0)]
)
def test_argument_validation(fast_client, query, by, limit):
    with pytest.raises(BadQuery):
        fast_client.lookup(query, by=by, limit=limit)  # type: ignore[arg-type]


def test_retries_on_503_then_succeeds(fast_client):
    with respx.mock(base_url="https://pubchem.ncbi.nlm.nih.gov/rest/pug") as router:
        route = router.get(path__regex=r"/compound/name/caffeine/property/.*")
        route.side_effect = [
            httpx.Response(503, json=fault("PUGREST.ServerBusy", "Too many requests")),
            httpx.Response(200, json=property_table(CAFFEINE_ROW)),
        ]
        assert fast_client.lookup("caffeine")[0].cid == 2519
        assert route.call_count == 2


def test_retries_exhausted_raises_unavailable(fast_client):
    with respx.mock(base_url="https://pubchem.ncbi.nlm.nih.gov/rest/pug") as router:
        route = router.get(path__regex=r"/compound/name/caffeine/property/.*")
        route.respond(503, json=fault("PUGREST.ServerBusy", "Too many requests"))
        with pytest.raises(PubChemUnavailable, match="busy"):
            fast_client.lookup("caffeine")
        assert route.call_count == fast_client.max_retries


def test_unexpected_status_and_non_json_body(fast_client):
    with respx.mock(base_url="https://pubchem.ncbi.nlm.nih.gov/rest/pug") as router:
        router.get(path__regex=r"/compound/name/x/property/.*").respond(
            500, text="<html>boom</html>"
        )
        with pytest.raises(PubChemError, match="HTTP 500"):
            fast_client.lookup("x")


def test_unexpected_payload_shape(fast_client):
    with respx.mock(base_url="https://pubchem.ncbi.nlm.nih.gov/rest/pug") as router:
        router.get(path__regex=r"/compound/name/x/property/.*").respond(
            200, json={"Something": "else"}
        )
        with pytest.raises(PubChemError, match="Unexpected"):
            fast_client.lookup("x")


def test_empty_property_table_is_not_found(fast_client):
    with respx.mock(base_url="https://pubchem.ncbi.nlm.nih.gov/rest/pug") as router:
        router.get(path__regex=r"/compound/name/x/property/.*").respond(
            200, json=property_table()
        )
        with pytest.raises(CompoundNotFound):
            fast_client.lookup("x")


def test_transport_error_is_wrapped(fast_client):
    with respx.mock(base_url="https://pubchem.ncbi.nlm.nih.gov/rest/pug") as router:
        router.get(path__regex=r"/compound/name/x/property/.*").mock(
            side_effect=httpx.ConnectError("no route")
        )
        with pytest.raises(PubChemError, match="failed"):
            fast_client.lookup("x")


def test_throttle_spaces_requests(monkeypatch):
    sleeps: list[float] = []
    monkeypatch.setattr("pubchem_finder.data.pug_rest.time.sleep", sleeps.append)
    clock = iter([0.0, 0.0, 0.05, 0.05, 10.0, 10.0])
    monkeypatch.setattr(
        "pubchem_finder.data.pug_rest.time.monotonic", lambda: next(clock)
    )
    with respx.mock(base_url="https://pubchem.ncbi.nlm.nih.gov/rest/pug") as router:
        router.get(path__regex=r"/compound/name/caffeine/property/.*").respond(
            200, json=property_table(CAFFEINE_ROW)
        )
        client = PubChemClient(min_interval=0.2)
        client.lookup("caffeine")
        client.lookup("caffeine")
        client.lookup("caffeine")
    # second call came 0.05s after the first -> sleep 0.15; third came 10s later -> no sleep
    assert sleeps == pytest.approx([0.15])


def test_throttle_serializes_threads():
    """Two threads can't both pass the throttle within min_interval: the MCP
    server calls the shared client from worker threads."""
    import threading

    client = PubChemClient(min_interval=0.2)
    client._throttle()  # prime _last_request so both threads have to wait
    barrier = threading.Barrier(2)
    passed: list[float] = []

    def go() -> None:
        barrier.wait()
        client._throttle()
        passed.append(time.monotonic())

    threads = [threading.Thread(target=go) for _ in range(2)]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
    client.close()
    assert len(passed) == 2
    # Without the lock both threads sleep at once and pass together. 30 ms of
    # slack covers Windows' timer granularity.
    assert abs(passed[1] - passed[0]) >= 0.2 - 0.03


def test_injected_httpx_client_is_not_closed():
    http = httpx.Client()
    client = PubChemClient(client=http)
    client.close()
    assert not http.is_closed
    http.close()


@pytest.mark.live
@pytest.mark.skipif(not os.environ.get("PUBCHEM_LIVE"), reason="set PUBCHEM_LIVE=1")
def test_live_lookup_caffeine():
    with PubChemClient() as client:
        rec = client.lookup("caffeine")[0]
    assert rec.cid == 2519
    assert rec.molecular_formula == "C8H10N4O2"


def test_unknown_properties_are_ignored_and_bad_values_rejected(fast_client):
    with respx.mock(base_url="https://pubchem.ncbi.nlm.nih.gov/rest/pug") as router:
        router.get(path__regex=r"/compound/name/x/property/.*").respond(
            200, json=property_table({**CAFFEINE_ROW, "NewProperty": 1})
        )
        assert fast_client.lookup("x")[0].cid == 2519
        router.get(path__regex=r"/compound/name/y/property/.*").respond(
            200, json=property_table({**CAFFEINE_ROW, "CID": "not-a-number"})
        )
        with pytest.raises(PubChemError, match="property values"):
            fast_client.lookup("y")
