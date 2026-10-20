"""Tests for the shared service layer."""

import pytest

from pubchem_finder import service
from pubchem_finder.core.exceptions import MoleculeValidationError


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
