"""Pytest configuration and fixtures."""

import pytest
from rdkit import Chem


@pytest.fixture
def ethanol_mol():
    """RDKit molecule object for ethanol (CCO)."""
    return Chem.MolFromSmiles("CCO")


@pytest.fixture
def ethanol_smiles():
    """SMILES string for ethanol."""
    return "CCO"


@pytest.fixture
def benzene_mol():
    """RDKit molecule object for benzene."""
    return Chem.MolFromSmiles("c1ccccc1")


@pytest.fixture
def benzene_smiles():
    """SMILES string for benzene."""
    return "c1ccccc1"


@pytest.fixture
def aspirin_smiles():
    """SMILES string for aspirin."""
    return "CC(=O)Oc1ccccc1C(=O)O"


@pytest.fixture
def aspirin_mol():
    """RDKit molecule object for aspirin."""
    return Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")


@pytest.fixture
def ethanol_with_salt_smiles():
    """SMILES string for ethanol with chloride salt."""
    return "CCO.Cl"


@pytest.fixture
def invalid_smiles():
    """Invalid SMILES string."""
    return "INVALID_SMILES_123"
