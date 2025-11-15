"""Tests for molecule standardizer."""

from rdkit import Chem

from pubchem_finder.chemistry import MoleculeStandardizer


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
    benzene_variants = ["c1ccccc1", "C1=CC=CC=C1", "c1cccc1c"]

    standardizer = MoleculeStandardizer()
    results = [standardizer.standardize_smiles(s) for s in benzene_variants]

    # All should give the same canonical SMILES
    assert len(set(results)) == 1
    assert results[0] is not None
