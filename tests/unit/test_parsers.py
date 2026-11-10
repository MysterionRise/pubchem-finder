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


def _assays(*rows):
    return {"Table": {"Columns": {"Column": fx.ASSAY_COLUMNS}, "Row": list(rows)}}


def test_most_potent_ranks_only_potency_measures():
    """The micromolar column also holds enzyme rates, plasma levels and solubilities.
    For quercetin, two Vmax values of a rate assay were listed as its 4th and 5th
    most potent activities."""
    summary = parsers.summarize_assays(
        _assays(
            fx.assay("200", "Active", "P1", "1", "0.006", "Vmax", "LOX-1 rate"),
            fx.assay("201", "Active", "P1", "1", "0.052", "Cp", "Plasma level"),
            fx.assay("202", "Active", "P1", "1", "2.3", "Solubility", "Solubility"),
            fx.assay("203", "Active", "P1", "1", "0.001", "", "Unnamed value"),
            fx.assay("204", "Active", "P1", "1", "0.5", "IC50", "X inhibition"),
            fx.assay("205", "Active", "P1", "1", "0.02", "ki", "X binding"),
        ),
        cid=1,
        top_n=10,
    )
    assert [(p.activity_name, p.value_um) for p in summary.most_potent] == [
        ("ki", 0.02),
        ("IC50", 0.5),
    ]
    assert summary.most_potent[0].target_gene_id == "1"
    assert summary.top_targets[0].active_assays == 6  # every active assay still counts


def test_targets_are_counted_per_gene_across_accessions():
    """PubChem lists one gene under a UniProt accession, NCBI protein IDs, or none.
    Each used to be its own target: ESR1 showed 11 active assays of quercetin's 14."""
    summary = parsers.summarize_assays(
        _assays(
            fx.assay("300", "Active", "NP_000116", "2099"),
            fx.assay("301", "Active", "P03372", "2099"),
            fx.assay("302", "Active", "", "2099"),
            fx.assay("302", "Active", "P03372", "2099"),  # one assay, two rows
            fx.assay("303", "Active", "P04637", "7157"),
            fx.assay("304", "Active", "P04637", "7157"),
            fx.assay("305", "Active", "XP_1", ""),  # no gene: the accession is the key
        ),
        cid=1,
    )
    assert [
        (t.target_gene_id, t.target_accession, t.active_assays)
        for t in summary.top_targets
    ] == [("2099", "P03372", 3), ("7157", "P04637", 2), (None, "XP_1", 1)]


def test_rows_without_a_gene_join_the_gene_their_accession_has_elsewhere():
    """ChEMBL's binding assays carry the accession and no gene, other depositors
    both. Imatinib showed ABL1 twice in its top five: 141 assays under gene 25 and
    49 under P00519 alone, where the table supports 190."""
    summary = parsers.summarize_assays(
        _assays(
            fx.assay("400", "Active", "P00519", "25"),
            fx.assay("401", "Active", "P00519", "", "0.01", "Ki", "ABL1 binding"),
            fx.assay("402", "Active", "P00519", ""),
            fx.assay("403", "Unspecified", "P16234", "5156"),  # the gene, not active
            fx.assay("404", "Active", "P16234", ""),
            # One accession under two genes: its gene-less row cannot be placed.
            fx.assay("405", "Active", "P00760", "615026"),
            fx.assay("406", "Inactive", "P00760", "780933"),
            fx.assay("407", "Active", "P00760", ""),
            fx.assay("408", "Active", "P00760", ""),
        ),
        cid=1,
        top_n=10,
    )
    assert [
        (t.target_gene_id, t.target_accession, t.active_assays)
        for t in summary.top_targets
    ] == [
        ("25", "P00519", 3),
        (None, "P00760", 2),
        ("5156", "P16234", 1),
        ("615026", "P00760", 1),
    ]
    assert [(p.aid, p.target_gene_id) for p in summary.most_potent] == [(401, "25")]


def test_target_accession_is_uniprots_even_when_only_other_rows_carry_it():
    """Quercetin's ESR1 was labelled AEP43755 (GenBank): PubChem lists P03372 for
    gene 2099 in the same table, on rows whose outcome is not active."""
    summary = parsers.summarize_assays(
        _assays(
            fx.assay("500", "Active", "AEP43755", "2099"),
            fx.assay("501", "Active", "AEP43755", "2099"),
            fx.assay("502", "Inactive", "P03372", "2099"),
            fx.assay("503", "Active", "NP_001", "7157"),  # no UniProt anywhere
        ),
        cid=1,
    )
    assert [(t.target_gene_id, t.target_accession) for t in summary.top_targets] == [
        ("2099", "P03372"),
        ("7157", "NP_001"),
    ]


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
