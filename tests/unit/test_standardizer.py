"""Tests for molecule standardizer."""

from pubchem_finder.chemistry import MoleculeStandardizer
from rdkit import Chem


def test_standardizer_removes_salts(ethanol_with_salt_smiles):
    """Test that standardization removes salt counterions."""
    standardizer = MoleculeStandardizer()
    result = standardizer.standardize_smiles(ethanol_with_salt_smiles)
    assert result == "CCO"


def test_standardizer_removes_hydrogens(ethanol_mol):
    """Test that standardization removes explicit hydrogens."""
    # Add explicit hydrogens
    mol_with_h = Chem.AddHs(ethanol_mol)
    assert mol_with_h.GetNumAtoms() > ethanol_mol.GetNumAtoms()

    standardizer = MoleculeStandardizer()
    std_mol = standardizer.standardize(mol_with_h)

    # Should have removed explicit H
    assert std_mol.GetNumAtoms() == ethanol_mol.GetNumAtoms()


def test_standardizer_largest_fragment():
    """Test that standardizer keeps largest fragment."""
    # Ethanol with sodium chloride salt
    smiles = "CCO.[Na+].[Cl-]"
    standardizer = MoleculeStandardizer()
    result = standardizer.standardize_smiles(smiles)
    assert result == "CCO"


def test_standardizer_neutralize_charges():
    """Test charge neutralization."""
    # Carboxylic acid anion
    smiles = "CC([O-])=O"
    standardizer = MoleculeStandardizer(neutralize=True)
    result = standardizer.standardize_smiles(smiles)

    # Should be neutralized to carboxylic acid
    mol = Chem.MolFromSmiles(result)
    total_charge = sum([atom.GetFormalCharge() for atom in mol.GetAtoms()])
    assert total_charge == 0


def test_standardizer_with_neutralize_false():
    """Test standardizer without neutralization."""
    smiles = "CC([O-])=O"
    standardizer = MoleculeStandardizer(neutralize=False)
    result = standardizer.standardize_smiles(smiles)

    mol = Chem.MolFromSmiles(result)
    total_charge = sum([atom.GetFormalCharge() for atom in mol.GetAtoms()])
    # Charge should be preserved
    assert total_charge == -1


def test_standardizer_with_invalid_smiles(invalid_smiles):
    """Test standardizer with invalid SMILES."""
    standardizer = MoleculeStandardizer()
    result = standardizer.standardize_smiles(invalid_smiles)
    assert result is None


def test_standardizer_with_none_mol():
    """Test standardizer with None molecule."""
    standardizer = MoleculeStandardizer()
    result = standardizer.standardize(None)
    assert result is None


def test_standardizer_canonicalization(benzene_smiles):
    """Test that different representations give same canonical form."""
    # Different valid SMILES representations of benzene
    benzene_variants = [
        "c1ccccc1",  # Aromatic notation
        "C1=CC=CC=C1",  # Kekulé form
        "c1cccc2c1ccc2",  # Naphthalene has multiple forms too
    ]

    standardizer = MoleculeStandardizer()

    # Test that aromatic and Kekulé benzene give same result
    benzene_aromatic = standardizer.standardize_smiles(benzene_variants[0])
    benzene_kekule = standardizer.standardize_smiles(benzene_variants[1])

    assert benzene_aromatic is not None
    assert benzene_kekule is not None
    assert benzene_aromatic == benzene_kekule
