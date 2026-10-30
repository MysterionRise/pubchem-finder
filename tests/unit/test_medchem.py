"""Tests for med-chem triage and scaffold grouping."""

import pytest
from rdkit import Chem

from pubchem_finder import medchem, service
from pubchem_finder.core.exceptions import MoleculeValidationError
from pubchem_finder.models import SimilarityHit

RHODANINE = "O=C1C(=Cc2ccccc2)SC(=S)N1"


def profile(smiles):
    return medchem.profile(Chem.MolFromSmiles(smiles), smiles)


def test_aspirin_profile():
    p = profile("CC(=O)Oc1ccccc1C(=O)O")
    assert 0.5 < p.qed < 0.6
    assert p.sa_score is not None
    assert 1 <= p.sa_score < 3
    assert p.lipinski_violations == []
    assert p.veber_ok is True
    assert p.pains_free is True
    assert [(a.catalog, a.name) for a in p.alerts] == [("Brenk", "phenol_ester")]
    assert p.murcko_scaffold == "c1ccccc1"


def test_pains_compound_flags_alerts_with_atoms():
    p = profile(RHODANINE)
    assert p.pains_free is False
    catalogs = {a.catalog for a in p.alerts}
    assert {"PAINS_A", "Brenk", "NIH"} <= catalogs
    pains = next(a for a in p.alerts if a.catalog == "PAINS_A")
    assert pains.atoms
    assert all(0 <= i < 14 for i in pains.atoms)
    assert pains.reference
    assert "Baell" in pains.reference


def test_rule_violations_and_acyclic_scaffold():
    p = profile("C" * 31 + "C(=O)O")
    assert "LogP > 5" in p.lipinski_violations
    assert p.veber_ok is False
    assert p.qed < 0.2
    assert p.murcko_scaffold == ""
    big = profile("O=C(O)" + "c1ccc(O)cc1" * 6)
    assert "molecular weight > 500" in big.lipinski_violations
    donors = profile("OCC(O)C(O)C(O)C(O)C(O)CO")
    assert "H-bond donors > 5" in donors.lipinski_violations
    acceptors = profile("O=C(OC)C(=O)OC(=O)OC(=O)OC(=O)OC(=O)OC(=O)OC")
    assert "H-bond acceptors > 10" in acceptors.lipinski_violations


def test_sa_score_missing_contrib(monkeypatch):
    monkeypatch.setattr(medchem, "_sascorer", lambda: None)
    assert medchem.sa_score(Chem.MolFromSmiles("CCO")) is None


def test_service_validates_smiles():
    with pytest.raises(MoleculeValidationError):
        service.medchem_profile("nope(")
    assert service.medchem_profile("CCO").canonical_smiles == "CCO"


def test_group_by_scaffold():
    hits = [
        SimilarityHit(cid=1, similarity=0.9, smiles="CC(=O)Oc1ccccc1C(=O)O"),
        SimilarityHit(cid=2, similarity=0.8, smiles="O=C(O)c1ccccc1O"),
        SimilarityHit(cid=3, similarity=0.85, smiles="O=C(O)c1ccc2ccccc2c1"),
        SimilarityHit(cid=4, similarity=0.5, smiles="CCO"),
        SimilarityHit(cid=5, similarity=0.4, smiles=None),
    ]
    groups = medchem.group_by_scaffold(hits)
    assert [g.scaffold for g in groups] == ["c1ccccc1", "c1ccc2ccccc2c1", "", "?"]
    assert groups[0].cids == [1, 2]
    assert groups[0].best_similarity == 0.9
