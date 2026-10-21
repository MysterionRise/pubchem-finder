"""Tests for the PubChem knowledge tools (client endpoints + service layer)."""

import os

import httpx
import pytest
import respx

from pubchem_finder import service
from pubchem_finder.core.exceptions import MoleculeValidationError
from pubchem_finder.data.pug_rest import (
    PROPERTIES,
    AnnotationNotFound,
    BadQuery,
    CompoundNotFound,
    PubChemClient,
    PubChemTimeout,
)
from tests import pubchem_fixtures as fx

ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"


# --------------------------------------------------------------------------- #
# Remote structure searches
# --------------------------------------------------------------------------- #


def test_similarity_search_request_and_ranking(pubchem_knowledge_api):
    result = service.pubchem_similarity_search(ASPIRIN, threshold=85, max_records=10)
    request = pubchem_knowledge_api.calls.last.request
    assert request.method == "POST"
    assert "/compound/fastsimilarity_2d/smiles/property/" + ",".join(PROPERTIES) in str(
        request.url
    )
    assert request.url.params["Threshold"] == "85"
    assert request.url.params["MaxRecords"] == "10"
    assert b"smiles=" in request.content

    assert result.search == "similarity"
    assert result.threshold == 85
    assert result.truncated is False
    cids = [h.cid for h in result.hits]
    assert cids[0] == 2244  # identical to the query: ECFP4 similarity 1.0
    assert result.hits[0].ecfp4_similarity == 1.0
    assert cids[-1] == 99999  # unparseable SMILES: no score, sorted last
    assert result.hits[-1].ecfp4_similarity is None
    scores = [h.ecfp4_similarity for h in result.hits[:-1]]
    assert scores == sorted(scores, reverse=True)
    assert "881-bit" in result.note
    assert result.hits[0].pubchem_url.endswith("/2244")


def test_truncated_when_max_records_reached(pubchem_knowledge_api):
    result = service.pubchem_similarity_search(ASPIRIN, max_records=4)
    assert result.truncated is True


def test_substructure_and_exact(pubchem_knowledge_api):
    sub = service.pubchem_substructure_search("c1ccccc1O", max_records=5)
    assert sub.search == "substructure"
    assert {h.cid for h in sub.hits} == {2244, 4133}
    exact = service.pubchem_exact_match(ASPIRIN, identity_type="same_connectivity")
    assert [h.cid for h in exact.hits] == [2244]
    assert pubchem_knowledge_api.calls.last.request.url.params["identity_type"] == (
        "same_connectivity"
    )
    assert "same_connectivity" in exact.note


def test_no_hits_is_an_empty_result(monkeypatch):
    monkeypatch.setattr("pubchem_finder.data.pug_rest.time.sleep", lambda s: None)
    with respx.mock() as router:
        router.post(url__regex=r".*fastsubstructure.*").respond(404, json=fx.NOT_FOUND)
        result = service.pubchem_substructure_search("C1CCCCC1")
    assert result.hits == []
    assert result.truncated is False


def test_invalid_query_rejected_before_any_request():
    with respx.mock() as router:
        with pytest.raises(MoleculeValidationError):
            service.pubchem_similarity_search("not smiles(")
        assert not router.calls


@pytest.mark.parametrize(
    ("call", "match"),
    [
        (lambda c: c.similarity_search(ASPIRIN, threshold=101), "threshold"),
        (lambda c: c.substructure_search(ASPIRIN, max_records=0), "max_records"),
        (lambda c: c.substructure_search(ASPIRIN, max_records=201), "max_records"),
        (lambda c: c.substructure_search("  "), "empty"),
        (lambda c: c.exact_match(ASPIRIN, identity_type="bogus"), "identity_type"),
        (lambda c: c.compound_json(0, "synonyms"), "positive"),
        (lambda c: c.view_data(0, "Toxicity"), "positive"),
        (lambda c: c.view_index(-1), "positive"),
    ],
)
def test_client_argument_validation(call, match):
    with PubChemClient(min_interval=0) as client, pytest.raises(BadQuery, match=match):
        call(client)


@pytest.mark.parametrize(
    "response",
    [
        httpx.Response(504, text="Gateway Timeout"),
        httpx.Response(503, json=fx.fault("PUGREST.Timeout", "Request timed out")),
        httpx.Response(500, json=fx.fault("PUGREST.Timeout", "Request timed out")),
    ],
)
def test_search_timeouts(response, monkeypatch):
    monkeypatch.setattr("pubchem_finder.data.pug_rest.time.sleep", lambda s: None)
    with respx.mock() as router:
        route = router.post(url__regex=r".*fastsubstructure.*").mock(
            return_value=response
        )
        with pytest.raises(PubChemTimeout):
            service.pubchem_substructure_search("c1ccccc1")
    assert route.call_count == 1  # a timeout is not retried


def test_transport_timeout(monkeypatch):
    monkeypatch.setattr("pubchem_finder.data.pug_rest.time.sleep", lambda s: None)
    with respx.mock() as router:
        router.post(url__regex=r".*fastsimilarity.*").mock(
            side_effect=httpx.ReadTimeout("slow")
        )
        with pytest.raises(PubChemTimeout, match="in time"):
            service.pubchem_similarity_search(ASPIRIN)


def test_invalid_json_body(monkeypatch):
    monkeypatch.setattr("pubchem_finder.data.pug_rest.time.sleep", lambda s: None)
    with respx.mock() as router:
        router.get(url__regex=r".*synonyms.*").respond(200, text="<html>")
        with pytest.raises(service.PubChemError, match="invalid JSON"):
            service.pubchem_synonyms(2244)


# --------------------------------------------------------------------------- #
# Per-compound information
# --------------------------------------------------------------------------- #


def test_synonyms(pubchem_knowledge_api):
    result = service.pubchem_synonyms(2244, limit=3)
    assert result.synonyms == ["aspirin", "ACETYLSALICYLIC ACID", "50-78-2"]
    assert result.cas_numbers == ["50-78-2"]
    assert result.total == 6
    with pytest.raises(BadQuery):
        service.pubchem_synonyms(2244, limit=0)
    with pytest.raises(CompoundNotFound):
        service.pubchem_synonyms(5)


def test_xrefs(pubchem_knowledge_api):
    result = service.pubchem_xrefs(2244, sample=1)
    assert result.counts == {"RegistryID": 3, "PubMedID": 2, "PatentID": 2}
    assert result.samples["RegistryID"] == ["50-78-2"]
    assert result.pubmed_urls == ["https://pubmed.ncbi.nlm.nih.gov/123/"]
    none = service.pubchem_xrefs(5)
    assert none.counts == {"RegistryID": 0, "PubMedID": 0, "PatentID": 0}
    with pytest.raises(BadQuery):
        service.pubchem_xrefs(2244, sample=0)


def test_description(pubchem_knowledge_api):
    result = service.pubchem_description(2244)
    assert result.title == "Aspirin"
    assert result.descriptions[0].source == "ChEBI"
    assert len(result.descriptions[1].text) <= 1500


def test_bioactivity(pubchem_knowledge_api):
    result = service.pubchem_bioactivity(2244)
    assert result.tested is True
    assert result.n_assays == 6
    assert result.top_targets[0].target_accession == "P23219"
    untested = service.pubchem_bioactivity(5)
    assert untested.tested is False
    assert untested.outcomes == {}


def test_annotations(pubchem_knowledge_api):
    result = service.pubchem_annotations(2244, "GHS Classification")
    request = pubchem_knowledge_api.calls.last.request
    assert request.url.params["heading"] == "GHS Classification"
    assert result.title == "Aspirin"
    assert result.n_total == 5
    assert result.truncated is False
    assert result.sources == ["European Chemicals Agency (ECHA)", "Some Depositor"]
    assert any(a.text == "Irritant" for a in result.annotations)

    trimmed = service.pubchem_annotations(2244, "GHS Classification", max_chars=500)
    assert trimmed.n_total == 5
    assert trimmed.truncated is True
    assert 1 <= len(trimmed.annotations) < 5

    numbers = service.pubchem_annotations(2244, "Experimental Properties")
    assert numbers.annotations[0].text == "135 °C"


def test_annotations_missing_section_lists_available(pubchem_knowledge_api):
    with pytest.raises(AnnotationNotFound, match="no 'Toxicity' section") as info:
        service.pubchem_annotations(2244, "Toxicity")
    assert "GHS Classification" in info.value.available
    assert "Safety and Hazards" in str(info.value)


def test_annotations_unknown_compound(pubchem_knowledge_api):
    with pytest.raises(CompoundNotFound, match="CID 7"):
        service.pubchem_annotations(7, "Toxicity")


def test_annotations_budget_validation():
    with pytest.raises(BadQuery):
        service.pubchem_annotations(2244, "Toxicity", max_chars=100)


@pytest.mark.parametrize(
    ("fn", "payload"),
    [
        (lambda: service.pubchem_synonyms(2244), {"x": 1}),
        (lambda: service.pubchem_xrefs(2244), {"x": 1}),
        (lambda: service.pubchem_description(2244), {"x": 1}),
        (lambda: service.pubchem_bioactivity(2244), {"x": 1}),
        (lambda: service.pubchem_annotations(2244, "Toxicity"), {"x": 1}),
    ],
)
def test_unexpected_payloads_become_pubchem_errors(fn, payload, monkeypatch):
    monkeypatch.setattr("pubchem_finder.data.pug_rest.time.sleep", lambda s: None)
    with respx.mock() as router:
        router.route().respond(200, json=payload)
        with pytest.raises(service.PubChemError, match="Unexpected"):
            fn()


# --------------------------------------------------------------------------- #
# Cache wiring in the shared client
# --------------------------------------------------------------------------- #


def test_shared_client_uses_cache_when_enabled(tmp_path, monkeypatch):
    monkeypatch.setenv("PUBCHEM_FINDER_HOME", str(tmp_path))
    monkeypatch.setenv("PUBCHEM_FINDER_CACHE", "1")
    monkeypatch.setattr(service, "_client", None)
    client = service.get_client()
    assert client.cache is not None
    assert client.cache.path == tmp_path / "cache" / "pubchem-http.sqlite"
    client.cache.put("k", b"{}")
    stats = service.cache_stats(tmp_path)
    assert stats is not None and stats.entries == 1
    assert service.cache_clear(tmp_path) == 1
    monkeypatch.setattr(service, "_client", None)


def test_shared_client_without_cache(monkeypatch, tmp_path):
    monkeypatch.setattr(service, "_client", None)
    assert service.get_client().cache is None  # PUBCHEM_FINDER_CACHE=0 in tests
    assert service.cache_stats(tmp_path / "none") is None
    assert service.cache_clear(tmp_path / "none") == 0
    monkeypatch.setattr(service, "_client", None)


def test_unwritable_cache_falls_back(monkeypatch, tmp_path):
    blocker = tmp_path / "file"
    blocker.write_text("not a directory")
    monkeypatch.setenv("PUBCHEM_FINDER_HOME", str(blocker))
    monkeypatch.setenv("PUBCHEM_FINDER_CACHE", "1")
    monkeypatch.setattr(service, "_client", None)
    assert service.get_client().cache is None
    monkeypatch.setattr(service, "_client", None)


# --------------------------------------------------------------------------- #
# Live checks (PUBCHEM_LIVE=1)
# --------------------------------------------------------------------------- #

live = pytest.mark.skipif(
    not os.environ.get("PUBCHEM_LIVE"), reason="set PUBCHEM_LIVE=1"
)


@pytest.mark.live
@live
def test_live_structure_searches():
    with PubChemClient() as client:
        similar = client.similarity_search(ASPIRIN, threshold=95, max_records=10)
        exact = client.exact_match(ASPIRIN)
        sub = client.substructure_search("CC(=O)Oc1ccccc1C(=O)O", max_records=5)
    assert 2244 in [r.cid for r in similar]
    assert [r.cid for r in exact][:1] == [2244]
    assert sub


@pytest.mark.live
@live
def test_live_compound_information(monkeypatch):
    monkeypatch.setattr(service, "_client", PubChemClient())
    assert "50-78-2" in service.pubchem_synonyms(2244).cas_numbers
    assert service.pubchem_xrefs(2244).counts["PubMedID"] > 0
    assert service.pubchem_description(2244).title == "Aspirin"
    bio = service.pubchem_bioactivity(2244)
    assert bio.tested and bio.n_assays > 0
    ghs = service.pubchem_annotations(2244, "GHS Classification")
    assert ghs.annotations and ghs.sources
    monkeypatch.setattr(service, "_client", None)
