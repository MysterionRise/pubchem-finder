"""Tests for the shared service layer."""

import pytest

from pubchem_finder import drawing, service
from pubchem_finder.core.exceptions import MoleculeValidationError


@pytest.mark.parametrize(
    "text",
    [
        "CCO CCN",
        "CCO\tCCN",
        "CCO\nCCN",
        "CCO ethanol",
        "CCO\N{GRINNING FACE}",
        "éCCO",
    ],
)
def test_text_after_a_space_or_non_ascii_is_refused_not_dropped(text):
    # RDKit reads 'SMILES name' and stops at non-ASCII, so these used to answer
    # as ethanol without a word about the rest of the input.
    with pytest.raises(MoleculeValidationError, match="one structure"):
        service.parse_smiles(text)
    with pytest.raises(MoleculeValidationError, match="one structure"):
        service.parse_query(text, "smarts")
    with pytest.raises(MoleculeValidationError, match="one structure"):
        drawing._parse([text])


def test_compare_refuses_second_smiles_instead_of_scoring_first_piece():
    with pytest.raises(MoleculeValidationError, match="one structure"):
        service.compare("CCO CCCCCCCC", "CCO")


def test_smarts_highlights_refuses_text_after_a_space():
    mol = drawing._parse(["CCO"])[0]
    with pytest.raises(MoleculeValidationError, match="one structure"):
        drawing.smarts_highlights([mol], "[#6] [#8]")


@pytest.mark.parametrize(
    "text",
    [
        "C[C@H](N)C(=O)O |&1:1|",  # a racemate, which came back as one enantiomer
        "C[C@H](N)C(=O)O |o1:1|",
        "*CC(*)c1ccccc1 |$star_e;;;star_e;;;;;;$,Sg:n:1,2,4,5,6,7,8,9::ht|",
    ],
)
def test_cxsmiles_is_refused_by_name_not_read_without_its_extension(text):
    # Issue #212: the stereo groups were parsed and then lost by every output.
    with pytest.raises(MoleculeValidationError, match="CXSMILES extension"):
        service.parse_smiles(text)
    with pytest.raises(MoleculeValidationError, match="CXSMILES extension"):
        drawing._parse([text])
    with pytest.raises(MoleculeValidationError, match="CXSMILES extension"):
        service.standardize(text)


PHENOL_AS_SMARTS = "[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#8]"


def test_smarts_atoms_are_refused_where_a_smiles_belongs():
    # Issue #212: read as SMILES these are atoms without hydrogens, and phenol
    # came back with a molecular weight of 88.07 (94.11) and no warning.
    for call in (
        lambda: service.properties(PHENOL_AS_SMARTS),
        lambda: service.medchem_profile(PHENOL_AS_SMARTS),
        lambda: service.compare(PHENOL_AS_SMARTS, "Oc1ccccc1"),
        lambda: service.parse_query(PHENOL_AS_SMARTS, "smiles"),
        lambda: drawing._parse(["CCO", PHENOL_AS_SMARTS]),
    ):
        with pytest.raises(MoleculeValidationError, match="atomic numbers") as err:
            call()
        assert 'query_type="smarts"' in str(err.value)
    with pytest.raises(MoleculeValidationError, match="SMILES at position 2"):
        drawing._parse(["CCO", PHENOL_AS_SMARTS])
    # RDKit writes a labelled atom with its isotope first, [13#6]: a SMARTS in
    # which every atom is labelled has no "[#" in it.
    for labelled in ("[13#6]1:[13#6]:[13#6]:[13#6]:[13#6]:[13#6]:1", "[18#8]"):
        with pytest.raises(MoleculeValidationError, match="atomic numbers"):
            service.properties(labelled)
    assert service.properties("[13cH]1[13cH][13cH][13cH][13cH][13cH]1").molecular_weight
    # As SMARTS it is a query, and still accepted as one.
    assert service.parse_query(PHENOL_AS_SMARTS, "smarts").GetNumAtoms() == 7
    # '#' as a triple bond next to a bracket atom is ordinary SMILES.
    assert service.properties("[C-]#[O+]").molecular_formula == "CO"
    assert service.properties("CC#N").molecular_formula == "C2H3N"


def test_properties_aspirin(aspirin_smiles):
    props = service.properties(aspirin_smiles)
    assert props.input_smiles == aspirin_smiles
    assert props.canonical_smiles == "CC(=O)Oc1ccccc1C(=O)O"
    assert props.molecular_weight == pytest.approx(180.16, abs=0.01)
    assert props.num_h_donors == 1
    assert props.num_aromatic_rings == 1
    assert props.lipinski_ok is True


def test_properties_invalid_smiles(invalid_smiles):
    with pytest.raises(MoleculeValidationError, match="Invalid SMILES"):
        service.properties(invalid_smiles)


@pytest.mark.parametrize("smiles", ["", "   "])
def test_empty_smiles_rejected(smiles):
    with pytest.raises(MoleculeValidationError):
        service.standardize(smiles)


def test_standardize_salt():
    result = service.standardize("CC(=O)[O-].[Na+]")
    assert result.smiles == "CC(=O)O"
    assert result.changed is True


def test_standardize_unchanged():
    result = service.standardize("CCO")
    assert result.smiles == "CCO"
    assert result.changed is False


def test_compare_identical_and_different(aspirin_smiles):
    same = service.compare(aspirin_smiles, "CC(=O)Oc1ccccc1C(=O)O")
    assert same.ecfp4_tanimoto == pytest.approx(1.0)
    assert same.maccs_tanimoto == pytest.approx(1.0)
    diff = service.compare(aspirin_smiles, "O=C(O)c1ccccc1O")
    assert 0.0 < diff.ecfp4_tanimoto < 1.0
    far = service.compare(aspirin_smiles, "C")
    assert far.ecfp4_tanimoto < diff.ecfp4_tanimoto


def test_compare_invalid_second_smiles(aspirin_smiles, invalid_smiles):
    with pytest.raises(MoleculeValidationError):
        service.compare(aspirin_smiles, invalid_smiles)


def test_lookup_uses_shared_client(pubchem_api, monkeypatch):
    monkeypatch.setattr("pubchem_finder.data.pug_rest.time.sleep", lambda s: None)
    result = service.lookup("caffeine")
    assert result.query == "caffeine"
    assert result.by == "name"
    assert result.compounds[0].cid == 2519
    assert service.get_client() is service.get_client()


# Oversized structures (hardening audit C1): RDKit's recursive SMILES writer overflows
# a thread's stack on long chains, which took the whole MCP server down.


def test_oversized_smiles_are_refused_before_rdkit_sees_them(monkeypatch):
    def parse(*_args, **_kwargs):
        raise AssertionError("parsed an oversized SMILES")

    monkeypatch.setattr(service.Chem, "MolFromSmiles", parse)
    with pytest.raises(MoleculeValidationError, match="80,000 characters"):
        service.properties("C" * 80_000)


@pytest.mark.parametrize(
    "call",
    [
        lambda s: service.properties(s),
        lambda s: service.standardize(s),
        lambda s: service.compare(s, "CCO"),
        lambda s: service.compare("CCO", s),
        lambda s: service.medchem_profile(s),
    ],
)
def test_structures_over_the_atom_limit_are_refused(call):
    with pytest.raises(
        MoleculeValidationError, match="1,001 atoms; the limit is 1,000"
    ):
        call("C" * 1_001)


def test_structures_at_the_limits_are_accepted():
    assert service.properties("C" * 1_000).num_heavy_atoms == 1_000


def test_oversized_smarts_queries_are_refused():
    with pytest.raises(MoleculeValidationError, match="SMARTS is 10,001 characters"):
        service.parse_query("C" * 10_001, "smarts")
    with pytest.raises(MoleculeValidationError, match="SMARTS has 1,001 atoms"):
        service.parse_query("[#6]" * 1_001, "smarts")


def test_get_client_keeps_one_client_per_data_directory(tmp_path, monkeypatch):
    """The client's cache follows its data directory, one client per directory."""
    monkeypatch.setattr(service, "_clients", {})
    monkeypatch.setenv("PUBCHEM_FINDER_CACHE", "1")
    home = tmp_path / "home"
    client = service.get_client(home)
    assert client.cache is not None
    assert client.cache.path == home.resolve() / "cache" / "pubchem-http.sqlite"
    assert service.get_client(home) is client
    assert service.get_client(tmp_path / "other") is not client


def test_compare_keeps_2048_bits_when_search_uses_1024(built_home, aspirin_smiles):
    comparison = service.compare(aspirin_smiles, "O=C(O)c1ccccc1O")
    results = service.similarity_search(aspirin_smiles, threshold=0.05, home=built_home)
    salicylic_acid = next(hit for hit in results.hits if hit.cid == 338)
    assert comparison.ecfp4_tanimoto == pytest.approx(0.4482758620689655)
    assert salicylic_acid.similarity == pytest.approx(0.4643)
