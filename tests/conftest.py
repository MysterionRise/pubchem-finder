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


# --------------------------------------------------------------------------- #
# Local index fixtures: small PubChem-style SDF chunks generated with RDKit.
# --------------------------------------------------------------------------- #

# (cid, smiles, iupac_name, smiles_tag) -- smiles_tag picks which PubChem tag
# carries the SMILES: "new" = PUBCHEM_SMILES, "old" = PUBCHEM_OPENEYE_ISO_SMILES,
# None = no SMILES tag (reader falls back to RDKit), "title" = no CID tag either,
# CID only in the title line.
CHUNK_A = [
    (2244, "CC(=O)Oc1ccccc1C(=O)O", "2-acetyloxybenzoic acid", "new"),
    (338, "O=C(O)c1ccccc1O", "2-hydroxybenzoic acid", "old"),
    (4133, "COC(=O)c1ccccc1O", "methyl 2-hydroxybenzoate", None),
    (243, "O=C(O)c1ccccc1", "benzoic acid", "title"),
    (702, "CCO", "ethanol", "new"),
]
CHUNK_B = [
    (2519, "Cn1cnc2c1c(=O)n(C)c(=O)n2C", "1,3,7-trimethylpurine-2,6-dione", "new"),
    (1983, "CC(=O)Nc1ccc(O)cc1", "N-(4-hydroxyphenyl)acetamide", "new"),
    (
        3672,
        "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
        "2-[4-(2-methylpropyl)phenyl]propanoic acid",
        "new",
    ),
    (
        156391,
        "COc1ccc2cc(ccc2c1)[C@H](C)C(=O)O",
        "(2S)-2-(6-methoxynaphthalen-2-yl)propanoic acid",
        "new",
    ),
    (4091, "CN(C)C(=N)N=C(N)N", "3-(diaminomethylidene)-1,1-dimethylguanidine", "new"),
    (241, "c1ccccc1", "benzene", "new"),
]
CHUNK_NAMES = (
    "Compound_000000001_000500000",
    "Compound_000500001_001000000",
    "Compound_001000001_001500000",
)
N_GOOD = len(CHUNK_A) + len(CHUNK_B)


def _good_block(cid, smiles, name, tag):
    from rdkit.Chem import Descriptors, rdMolDescriptors

    mol = Chem.MolFromSmiles(smiles)
    mol.SetProp("_Name", str(cid) if tag == "title" else "")
    if tag != "title":
        mol.SetProp("PUBCHEM_COMPOUND_CID", str(cid))
    if tag == "new":
        mol.SetProp("PUBCHEM_SMILES", smiles)
    elif tag == "old":
        mol.SetProp("PUBCHEM_OPENEYE_ISO_SMILES", smiles)
    mol.SetProp("PUBCHEM_IUPAC_NAME", name)
    mol.SetProp("PUBCHEM_IUPAC_INCHIKEY", Chem.MolToInchiKey(mol))
    mol.SetProp("PUBCHEM_MOLECULAR_FORMULA", rdMolDescriptors.CalcMolFormula(mol))
    mol.SetProp("PUBCHEM_MOLECULAR_WEIGHT", f"{Descriptors.MolWt(mol):.2f}")
    return mol


def _bad_valence_block(cid):
    mol = Chem.MolFromSmiles("C(C)(C)(C)(C)C", sanitize=False)
    mol.UpdatePropertyCache(strict=False)
    return Chem.MolToMolBlock(mol) + f">  <PUBCHEM_COMPOUND_CID>\n{cid}\n\n$$$$\n"


def write_sdf(path, records, extra_blocks=(), gz=True):
    """Write PubChem-style records (see CHUNK_A) to an SDF or SDF.gz file."""
    import gzip
    import io

    buf = io.StringIO()
    writer = Chem.SDWriter(buf)
    for rec in records:
        writer.write(_good_block(*rec))
    writer.close()
    text = buf.getvalue() + "".join(extra_blocks)
    if gz:
        with gzip.open(path, "wt") as fh:
            fh.write(text)
    else:
        path.write_text(text)
    return path


def no_cid_block():
    mol = Chem.MolFromSmiles("CCN")
    mol.SetProp("_Name", "not-a-cid")
    import io

    buf = io.StringIO()
    writer = Chem.SDWriter(buf)
    writer.write(mol)
    writer.close()
    return buf.getvalue()


@pytest.fixture(scope="session")
def sample_sdfs(tmp_path_factory):
    """Three gzipped chunks: A (5 good, 2 bad), B (6 good), C (only bad)."""
    d = tmp_path_factory.mktemp("sdf")
    a = write_sdf(
        d / f"{CHUNK_NAMES[0]}.sdf.gz",
        CHUNK_A,
        extra_blocks=(_bad_valence_block(999999), no_cid_block()),
    )
    b = write_sdf(d / f"{CHUNK_NAMES[1]}.sdf.gz", CHUNK_B)
    c = write_sdf(
        d / f"{CHUNK_NAMES[2]}.sdf.gz", [], extra_blocks=(_bad_valence_block(1),)
    )
    return a, b, c


@pytest.fixture(scope="session")
def built_home(tmp_path_factory, sample_sdfs):
    """A data directory with a finalized index over chunks A, B and C."""
    from pubchem_finder.index.builder import SdfSource, build_chunks, finalize
    from pubchem_finder.index.layout import IndexPaths, chunk_name

    home = tmp_path_factory.mktemp("home")
    paths = IndexPaths(home)
    sources = [SdfSource(name=chunk_name(p), local_path=p) for p in sample_sdfs]
    build_chunks(paths, sources, workers=1)
    finalize(paths)
    return home


@pytest.fixture(scope="session")
def substructure_home(tmp_path_factory, sample_sdfs):
    """Like built_home, with substructure search enabled."""
    from pubchem_finder.index.builder import SdfSource, build_chunks, finalize
    from pubchem_finder.index.layout import IndexPaths, chunk_name

    home = tmp_path_factory.mktemp("subhome")
    paths = IndexPaths(home)
    sources = [SdfSource(name=chunk_name(p), local_path=p) for p in sample_sdfs]
    build_chunks(paths, sources, workers=1)
    finalize(paths, substructure=True)
    return home


def pytest_configure(config):
    """Keep tests away from the user's real data directory and response cache."""
    import os
    import tempfile

    os.environ["PUBCHEM_FINDER_HOME"] = tempfile.mkdtemp(prefix="pubchem-finder-tests-")
    os.environ["PUBCHEM_FINDER_CACHE"] = "0"


@pytest.fixture
def pubchem_knowledge_api(monkeypatch):
    """respx mock of the Phase 4a endpoints, with retry/throttle sleeps disabled."""
    import respx

    from tests.pubchem_fixtures import mock_pubchem

    monkeypatch.setattr("pubchem_finder.data.pug_rest.time.sleep", lambda s: None)
    with respx.mock(assert_all_called=False) as router:
        yield mock_pubchem(router)
