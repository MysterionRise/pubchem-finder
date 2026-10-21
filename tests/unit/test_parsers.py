"""Tests for the pure PubChem payload parsers."""

import pytest

from pubchem_finder.data import parsers
from pubchem_finder.models import Annotation
from tests import pubchem_fixtures as fx


@pytest.mark.parametrize(
    ("text", "valid"),
    [
        ("50-78-2", True),  # aspirin
        ("7732-18-5", True),  # water
        ("58-08-2", True),  # caffeine
        ("50-78-3", False),  # wrong check digit
        ("2024-01-15", False),  # a date
        ("1-23-4", False),  # too short
        ("aspirin", False),
    ],
)
def test_is_cas_number(text, valid):
    assert parsers.is_cas_number(text) is valid


def test_parse_synonyms_dedupes_case_insensitively():
    names = parsers.parse_synonyms(fx.SYNONYMS)
    assert names[0] == "aspirin"
    assert "Aspirin" not in names
    assert len(names) == 6


def test_parse_xrefs_dedupes_and_keeps_order():
    refs = parsers.parse_xrefs(
        fx.XREFS, ("RegistryID", "PubMedID", "PatentID", "Missing")
    )
    assert refs["RegistryID"] == ["50-78-2", "CHEMBL25", "DB00945"]
    assert refs["PubMedID"] == ["123", "456"]
    assert refs["Missing"] == []


def test_parse_descriptions_title_sources_and_truncation():
    title, texts = parsers.parse_descriptions(fx.DESCRIPTION, max_chars=100)
    assert title == "Aspirin"
    assert texts[0].source == "ChEBI"
    assert texts[0].url.startswith("https://www.ebi.ac.uk/")
    assert len(texts[1].text) == 100
    assert texts[1].text.endswith("…")
    _, one = parsers.parse_descriptions(fx.DESCRIPTION, max_items=1)
    assert len(one) == 1


@pytest.mark.parametrize(
    "bad", [{}, {"InformationList": {"Information": "nope"}}, None]
)
def test_information_payload_errors(bad):
    with pytest.raises(parsers.PayloadError):
        parsers.parse_synonyms(bad)


def test_summarize_assays():
    summary = parsers.summarize_assays(fx.ASSAY_SUMMARY, cid=2244)
    assert summary.tested is True
    assert summary.n_assays == 6  # the row with a non-numeric AID is ignored
    assert summary.outcomes == {"Active": 5, "Inactive": 1, "Inconclusive": 1}
    top = summary.top_targets[0]
    assert (top.target_accession, top.target_gene_id, top.active_assays) == (
        "P23219",
        "5742",
        2,
    )
    assert summary.top_targets[1].target_accession == "P35354"
    assert [p.value_um for p in summary.most_potent] == [0.5, 0.7, 1.67, 278.0]
    assert summary.most_potent[0].aid == 102
    assert summary.pubchem_url.endswith("/2244#section=BioAssay-Results")


def test_summarize_assays_empty_and_malformed():
    empty = {"Table": {"Columns": {"Column": fx.ASSAY_COLUMNS}, "Row": []}}
    summary = parsers.summarize_assays(empty, cid=1)
    assert summary.tested is False
    assert summary.n_assays == 0
    with pytest.raises(parsers.PayloadError, match="AID/outcome"):
        parsers.summarize_assays({"Table": {"Columns": {"Column": ["X"]}}}, cid=1)
    with pytest.raises(parsers.PayloadError):
        parsers.summarize_assays({"Nope": 1}, cid=1)


def test_parse_annotations_ghs():
    title, annotations = parsers.parse_annotations(fx.GHS_VIEW)
    assert title == "Aspirin"
    by_name = {(a.name, a.source): a for a in annotations}
    pictogram = by_name[("Pictogram(s)", "European Chemicals Agency (ECHA)")]
    assert pictogram.text == "Irritant"  # meaning comes from the icon markup
    assert (
        pictogram.section
        == "Safety and Hazards > Hazards Identification > GHS Classification"
    )
    assert pictogram.license == "https://echa.europa.eu/web/guest/legal-notice"
    hazards = by_name[("GHS Hazard Statements", "European Chemicals Agency (ECHA)")]
    assert (
        hazards.text
        == "H302 (100%): Harmful if swallowed; H315: Causes skin irritation"
    )
    assert by_name[("Signal", "Some Depositor")].text == "Danger"
    assert all(a.name != "Empty" for a in annotations)  # empty values are dropped


def test_parse_annotations_numbers_and_units():
    _, annotations = parsers.parse_annotations(fx.MELTING_POINT_VIEW)
    assert [a.text for a in annotations] == ["135 °C", "135.5, 136"]
    assert annotations[0].section == "Experimental Properties > Melting Point"
    assert annotations[0].source == "HSDB"


def test_parse_headings():
    assert parsers.parse_headings(fx.VIEW_INDEX) == [
        "Names and Identifiers",
        "Computed Descriptors",
        "Safety and Hazards",
        "Hazards Identification",
        "GHS Classification",
    ]
    with pytest.raises(parsers.PayloadError):
        parsers.parse_headings({})
    with pytest.raises(parsers.PayloadError):
        parsers.parse_annotations({})


def test_trim_annotations():
    items = [Annotation(section="S", text="x" * 100) for _ in range(5)]
    kept, truncated = parsers.trim_annotations(items, max_chars=300)
    assert len(kept) == 2
    assert truncated is True
    kept, truncated = parsers.trim_annotations(items, max_chars=10_000)
    assert len(kept) == 5
    assert truncated is False
    huge = [Annotation(section="S", text="y" * 5000)]
    kept, truncated = parsers.trim_annotations(huge, max_chars=1000)
    assert truncated is True
    assert len(kept[0].text) < 1000
    assert kept[0].text.endswith("…")
