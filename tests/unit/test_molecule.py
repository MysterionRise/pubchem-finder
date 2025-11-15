"""Tests for core Molecule data structure."""

import pytest

from pubchem_finder.core import Molecule, MoleculeValidationError


def test_molecule_creation_with_smiles():
    """Test creating a molecule with SMILES."""
    mol = Molecule(smiles="CCO", canonical_smiles="CCO")
    assert mol.smiles == "CCO"
    assert mol.canonical_smiles == "CCO"


def test_molecule_requires_smiles():
    """Test that molecule requires either smiles or canonical_smiles."""
    with pytest.raises(MoleculeValidationError):
        Molecule()


def test_molecule_with_only_canonical_smiles():
    """Test creating molecule with only canonical SMILES."""
    mol = Molecule(canonical_smiles="CCO")
    assert mol.canonical_smiles == "CCO"


def test_molecule_add_fingerprint():
    """Test adding fingerprint to molecule."""
    mol = Molecule(smiles="CCO")
    mol.add_fingerprint("morgan", [1, 0, 1, 0])
    assert mol.has_fingerprint("morgan")
    assert mol.get_fingerprint("morgan") == [1, 0, 1, 0]


def test_molecule_get_nonexistent_fingerprint():
    """Test getting fingerprint that doesn't exist."""
    mol = Molecule(smiles="CCO")
    assert mol.get_fingerprint("morgan") is None
    assert not mol.has_fingerprint("morgan")


def test_molecule_with_properties():
    """Test molecule with molecular properties."""
    mol = Molecule(
        smiles="CCO",
        canonical_smiles="CCO",
        molecular_weight=46.07,
        logp=-0.07,
        num_h_donors=1,
        num_h_acceptors=1,
    )
    assert mol.molecular_weight == 46.07
    assert mol.logp == -0.07
    assert mol.num_h_donors == 1
    assert mol.num_h_acceptors == 1


def test_molecule_repr():
    """Test molecule string representation."""
    mol = Molecule(smiles="CCO", canonical_smiles="CCO", molecular_weight=46.07)
    repr_str = repr(mol)
    assert "CCO" in repr_str
    assert "46.07" in repr_str


def test_molecule_with_id():
    """Test molecule with ID."""
    mol = Molecule(id="PubChem_12345", smiles="CCO")
    assert mol.id == "PubChem_12345"


def test_molecule_created_at():
    """Test that created_at is set automatically."""
    mol = Molecule(smiles="CCO")
    assert mol.created_at is not None
