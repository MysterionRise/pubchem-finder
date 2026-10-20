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


@pytest.fixture
def anyio_backend():
    """Run async tests on asyncio only (the MCP server is asyncio-based)."""
    return "asyncio"


CAFFEINE_ROW = {
    "CID": 2519,
    "Title": "Caffeine",
    "IUPACName": "1,3,7-trimethylpurine-2,6-dione",
    "SMILES": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "ConnectivitySMILES": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "InChIKey": "RYYVLZVUVIJVGH-UHFFFAOYSA-N",
    "MolecularFormula": "C8H10N4O2",
    "MolecularWeight": "194.19",  # PUG REST serializes MolecularWeight as a string
    "XLogP": -0.1,
    "HBondDonorCount": 0,
    "HBondAcceptorCount": 3,
    "TPSA": 61.8,
}

ASPIRIN_ROW = {
    "CID": 2244,
    "Title": "Aspirin",
    "IUPACName": "2-acetyloxybenzoic acid",
    "SMILES": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "ConnectivitySMILES": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "InChIKey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
    "MolecularFormula": "C9H8O4",
    "MolecularWeight": "180.16",
    "XLogP": 1.2,
    "HBondDonorCount": 1,
    "HBondAcceptorCount": 4,
    "TPSA": 63.6,
}


def property_table(*rows):
    return {"PropertyTable": {"Properties": list(rows)}}


def fault(code: str, message: str, *details: str):
    return {"Fault": {"Code": code, "Message": message, "Details": list(details)}}


@pytest.fixture
def pubchem_api():
    """A respx router mocking PUG REST, with caffeine/aspirin canned answers."""
    import respx

    base = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
    with respx.mock(base_url=base, assert_all_called=False) as router:
        router.get(path__regex=r"/compound/name/caffeine/property/.*").respond(
            200, json=property_table(CAFFEINE_ROW)
        )
        router.get(path__regex=r"/compound/cid/2244/property/.*").respond(
            200, json=property_table(ASPIRIN_ROW)
        )
        router.get(path__regex=r"/compound/cid/2244,2519/property/.*").respond(
            200, json=property_table(ASPIRIN_ROW, CAFFEINE_ROW)
        )
        router.get(
            path__regex=r"/compound/inchikey/BSYNRYMUTXBXSQ-UHFFFAOYSA-N/property/.*"
        ).respond(200, json=property_table(ASPIRIN_ROW))
        router.post(path__regex=r"/compound/smiles/property/.*").respond(
            200, json=property_table(ASPIRIN_ROW)
        )
        router.get(path__regex=r"/compound/name/notacompound/property/.*").respond(
            404,
            json=fault(
                "PUGREST.NotFound",
                "No CID found",
                "No CID found that matches the given name",
            ),
        )
        router.get(path__regex=r"/compound/cid/abc/property/.*").respond(
            400, json=fault("PUGREST.BadRequest", "Invalid CID")
        )
        yield router


@pytest.fixture
def fast_client():
    """A PubChemClient with no throttling or retry sleeps."""
    from unittest.mock import patch

    from pubchem_finder.data.pug_rest import PubChemClient

    with (
        patch("pubchem_finder.data.pug_rest.time.sleep"),
        PubChemClient(min_interval=0) as client,
    ):
        yield client
