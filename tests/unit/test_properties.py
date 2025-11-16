"""Tests for molecular property calculation."""

from pubchem_finder.chemistry import calculate_properties
from pubchem_finder.chemistry.properties import is_druglike
from rdkit import Chem


def test_calculate_properties_ethanol(ethanol_mol):
    """Test property calculation for ethanol."""
    props = calculate_properties(ethanol_mol)

    assert "molecular_weight" in props
    assert "logp" in props
    assert "num_h_donors" in props
    assert "num_h_acceptors" in props
    assert "tpsa" in props

    # Ethanol (CCO) has MW ~46, 1 H-donor, 1 H-acceptor
    assert abs(props["molecular_weight"] - 46.07) < 0.1
    assert props["num_h_donors"] == 1
    assert props["num_h_acceptors"] == 1


def test_calculate_properties_benzene(benzene_mol):
    """Test property calculation for benzene."""
    props = calculate_properties(benzene_mol)

    # Benzene (C6H6) has MW ~78, no H-donors/acceptors
    assert abs(props["molecular_weight"] - 78.11) < 0.1
    assert props["num_h_donors"] == 0
    assert props["num_h_acceptors"] == 0
    assert props["num_aromatic_rings"] == 1


def test_calculate_properties_aspirin(aspirin_mol):
    """Test property calculation for aspirin."""
    props = calculate_properties(aspirin_mol)

    # Aspirin has 2 H-acceptors (from carboxylic acid)
    assert props["num_h_donors"] >= 1
    assert props["num_h_acceptors"] >= 2
    assert props["molecular_weight"] > 100


def test_is_druglike_ethanol(ethanol_mol):
    """Test drug-likeness check for ethanol (should be drug-like)."""
    assert is_druglike(ethanol_mol)


def test_is_druglike_large_molecule():
    """Test drug-likeness for large molecule (should fail)."""
    # Create a large molecule that violates Lipinski's rules
    # Very large molecule with high MW
    large_smiles = "C" * 100  # Long alkane chain
    mol = Chem.MolFromSmiles(large_smiles)

    assert not is_druglike(mol)


def test_properties_types(ethanol_mol):
    """Test that all properties are correct types."""
    props = calculate_properties(ethanol_mol)

    assert isinstance(props["molecular_weight"], float)
    assert isinstance(props["logp"], float)
    assert isinstance(props["num_h_donors"], int)
    assert isinstance(props["num_h_acceptors"], int)
    assert isinstance(props["tpsa"], float)
    assert isinstance(props["num_rotatable_bonds"], int)
