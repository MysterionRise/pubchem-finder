"""Tests for the streaming SDF reader."""

import pytest
from rdkit import Chem

from pubchem_finder.data.sdf import SdfReader
from tests.conftest import CHUNK_A, write_sdf


def test_reads_good_records_and_counts_failures(sample_sdfs):
    reader = SdfReader(sample_sdfs[0])
    records = {r.cid: r for r in reader}
    assert set(records) == {cid for cid, *_ in CHUNK_A}
    assert reader.n_ok == len(CHUNK_A)
    assert reader.n_failed == 2  # over-valent carbon, and a record with no CID


def test_tag_fallbacks(sample_sdfs):
    records = {r.cid: r for r in SdfReader(sample_sdfs[0])}
    assert records[2244].smiles == "CC(=O)Oc1ccccc1C(=O)O"  # PUBCHEM_SMILES
    assert records[338].smiles == "O=C(O)c1ccccc1O"  # old OpenEye tag
    assert records[4133].smiles == Chem.MolToSmiles(
        Chem.MolFromSmiles("COC(=O)c1ccccc1O")
    )
    assert 243 in records  # CID taken from the title line
    aspirin = records[2244]
    assert aspirin.iupac_name == "2-acetyloxybenzoic acid"
    assert aspirin.inchikey == "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"
    assert aspirin.molecular_formula == "C9H8O4"
    assert aspirin.molecular_weight == pytest.approx(180.16)


def test_uncompressed_file_and_bad_weight(tmp_path):
    path = write_sdf(tmp_path / "plain.sdf", [(702, "CCO", "ethanol", "new")], gz=False)
    text = path.read_text().replace("46.07", "n/a")
    path.write_text(text)
    (record,) = list(SdfReader(path))
    assert record.cid == 702
    assert record.molecular_weight is None


def test_empty_file(tmp_path):
    path = tmp_path / "empty.sdf"
    path.write_text("")
    reader = SdfReader(path)
    assert list(reader) == []
    assert (reader.n_ok, reader.n_failed) == (0, 0)
