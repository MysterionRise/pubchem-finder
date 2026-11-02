"""The shipped demo set, and scripts/build_demo_set.py that regenerates it."""

import csv
import importlib.util
import json
import sys
from pathlib import Path

import respx
from rdkit import Chem

from pubchem_finder.index.builder import METADATA_SCHEMA

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "src" / "pubchem_finder" / "index" / "demo_data"
TSV = DATA / "fda_approved_drugs.tsv"
PUG = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
VIEW = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view"


def load_script():
    spec = importlib.util.spec_from_file_location(
        "build_demo_set", ROOT / "scripts" / "build_demo_set.py"
    )
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def read_rows():
    with TSV.open(newline="", encoding="utf-8") as fh:
        return list(csv.reader(fh, delimiter="\t"))


def test_shipped_tsv_has_the_index_columns_one_row_per_cid():
    header, *rows = read_rows()
    assert header == METADATA_SCHEMA.names
    cids = [int(r[0]) for r in rows]
    assert cids == sorted(set(cids))
    assert all(len(r) == len(header) for r in rows)
    meta = json.loads((DATA / "fda_approved_drugs.json").read_text(encoding="utf-8"))
    assert meta["n_compounds"] == len(rows) > 1500
    assert "public domain" in meta["license"]


def test_shipped_tsv_is_lf_only():
    """Bytes matter: the demo index is keyed by a hash of the file."""
    data = TSV.read_bytes()
    assert b"\r" not in data
    assert data.endswith(b"\n")


def test_every_shipped_structure_parses_as_one_component():
    _, *rows = read_rows()
    for cid, smiles, *_ in rows:
        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None, cid
        assert "." not in smiles, cid


def test_shipped_set_has_the_drugs_the_docs_use():
    _, *rows = read_rows()
    cids = {int(r[0]) for r in rows}
    # aspirin, acetaminophen, ibuprofen, caffeine, metformin (not its hydrochloride)
    assert {2244, 1983, 3672, 2519, 4091} <= cids
    assert 14219 not in cids


def annotation(name, cids, source="Drugs@FDA"):
    return {
        "SourceName": source,
        "Name": name,
        "URL": "https://www.accessdata.fda.gov/scripts/cder/daf/",
        "LicenseNote": "in the public domain",
        "LicenseURL": "https://www.fda.gov/about-fda/about-website/website-policies",
        "LinkedRecords": {"CID": cids},
    }


def row(cid, smiles, name="x"):
    return {
        "CID": cid,
        "SMILES": smiles,
        "InChIKey": f"KEY{cid}",
        "IUPACName": name,
        "MolecularFormula": "C",
        "MolecularWeight": "1.5",
    }


def test_generator_writes_parents_and_provenance(tmp_path, monkeypatch):
    script = load_script()
    monkeypatch.setattr("pubchem_finder.data.pug_rest.time.sleep", lambda s: None)
    pages = {
        "1": [
            annotation("ASPIRIN", [2244]),
            annotation("METFORMIN HYDROCHLORIDE", [14219]),
            annotation("A REVIEW", [999], source="Someone's annotations"),
        ],
        "2": [
            annotation("LEVOCARNITINE", [10917]),
            annotation("A MIXTURE", [555]),
            annotation("NO STRUCTURE", []),
        ],
    }
    with respx.mock(assert_all_called=True) as router:
        router.get(f"{VIEW}/annotations/heading/JSON").mock(
            side_effect=lambda request: respx.MockResponse(
                200,
                json={
                    "Annotations": {
                        "Annotation": pages[request.url.params["page"]],
                        "TotalPages": 2,
                    }
                },
            )
        )
        parents = router.post(f"{PUG}/compound/cid/cids/JSON").respond(
            200, json={"IdentifierList": {"CID": [2244, 4091]}}
        )
        props = router.post(path__regex=r"/rest/pug/compound/cid/property/.*").respond(
            200,
            json={
                "PropertyTable": {
                    "Properties": [
                        row(2244, "CC(=O)Oc1ccccc1C(=O)O", "2-acetyloxybenzoic acid"),
                        row(4091, "CN(C)C(=N)N=C(N)N"),
                        row(14219, "CN(C)C(=N)N=C(N)N.Cl"),
                        row(10917, "C[N+](C)(C)C[C@@H](CC(=O)[O-])O", ""),
                        row(555, "CCO.CCN"),
                    ]
                }
            },
        )
        with script.PubChemClient(min_interval=0) as client:
            n_rows = script.build(client, tmp_path, "2026-11-01")

    assert parents.calls.last.request.url.params["cids_type"] == "parent"
    assert parents.calls.last.request.content == b"cid=555%2C2244%2C10917%2C14219"
    assert props.calls.last.request.content == (
        b"cid=555%2C2244%2C4091%2C10917%2C14219"
    )
    assert n_rows == 3
    assert (tmp_path / "fda_approved_drugs.tsv").read_text() == (
        "cid\tsmiles\tinchikey\tiupac_name\tmolecular_formula\tmolecular_weight\n"
        "2244\tCC(=O)Oc1ccccc1C(=O)O\tKEY2244\t2-acetyloxybenzoic acid\tC\t1.5\n"
        "4091\tCN(C)C(=N)N=C(N)N\tKEY4091\tx\tC\t1.5\n"
        "10917\tC[N+](C)(C)C[C@@H](CC(=O)[O-])O\tKEY10917\t\tC\t1.5\n"
    )
    meta = json.loads((tmp_path / "fda_approved_drugs.json").read_text())
    assert meta["n_compounds"] == 3
    assert meta["retrieved"] == "2026-11-01"
    assert meta["license"] == "in the public domain"
    assert "the 4 compounds PubChem links to Drugs@FDA" in meta["description"]
