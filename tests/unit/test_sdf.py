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


# FPSim2 keeps compound IDs as unsigned 32-bit numbers. A larger ID came back from a
# search as that ID modulo 2**32: another compound's number, at similarity 1.0.


def test_ids_above_32_bits_are_skipped_not_wrapped(tmp_path):
    path = write_sdf(
        tmp_path / "ids.sdf.gz",
        [
            (2244, "CC(=O)Oc1ccccc1C(=O)O", "aspirin", "new"),
            (2**32 + 702, "CCO", "wraps to 702", "new"),
            (2**32 - 1, "CCN", "the largest that fits", "new"),
            (0, "CCC", "zero is no ID", "new"),
            (241, "c1ccccc1", "benzene", "title"),
        ],
    )
    reader = SdfReader(path)
    assert [r.cid for r in reader] == [2244, 2**32 - 1, 241]
    assert reader.n_failed == 2
    assert reader.skipped_id == str(2**32 + 702)


@pytest.mark.parametrize("title", ["\u00b2", "\u0663", "9" * 5000, "12a"])
def test_ids_that_are_not_plain_numbers_are_failed_records(tmp_path, title):
    """str.isdigit() accepts "\u00b2" (int() does not) and "\u0663" (int() reads 3)."""
    mol = Chem.MolFromSmiles("CCO")
    mol.SetProp("_Name", title)
    path = tmp_path / "odd.sdf"
    with Chem.SDWriter(str(path)) as writer:
        writer.write(mol)
    reader = SdfReader(path)
    assert list(reader) == []
    assert reader.n_failed == 1
    # The first unusable ID is kept for the build's message, shortened.
    assert reader.skipped_id == (title if len(title) <= 24 else "9" * 21 + "...")


@pytest.mark.parametrize(
    ("text", "cid"),
    [("000000002244", 2244), ("0" * 40 + "7", 7), ("04294967295", 2**32 - 1)],
)
def test_zero_padded_ids_are_read_as_their_number(tmp_path, text, cid):
    """The length guard used to run before the zeros came off, so a registry that
    pads its numbers to 11 or more digits lost every record."""
    mol = Chem.MolFromSmiles("CCO")
    mol.SetProp("_Name", text)
    path = tmp_path / "padded.sdf"
    with Chem.SDWriter(str(path)) as writer:
        writer.write(mol)
    reader = SdfReader(path)
    assert [r.cid for r in reader] == [cid]
    assert reader.skipped_id is None


@pytest.mark.parametrize("text", ["0", "0000", "04294967296"])
def test_zero_and_padded_out_of_range_ids_are_skipped(tmp_path, text):
    mol = Chem.MolFromSmiles("CCO")
    mol.SetProp("_Name", text)
    path = tmp_path / "zeros.sdf"
    with Chem.SDWriter(str(path)) as writer:
        writer.write(mol)
    reader = SdfReader(path)
    assert list(reader) == []
    assert reader.skipped_id == text


def _tagged(path, tag, title, smiles="CCO"):
    mol = Chem.MolFromSmiles(smiles)
    mol.SetProp("_Name", title)
    mol.SetProp("PUBCHEM_COMPOUND_CID", tag)
    with Chem.SDWriter(str(path)) as writer:
        writer.write(mol)
    return path


def test_out_of_range_tag_is_not_replaced_by_a_numeric_title(tmp_path):
    """A registry number above 32 bits in the tag and a row counter as the title:
    the record was indexed as compound 7, with nothing said."""
    reader = SdfReader(_tagged(tmp_path / "t.sdf", str(2**32 + 702), "7"))
    assert list(reader) == []
    assert reader.n_failed == 1
    assert reader.skipped_id == str(2**32 + 702)


def test_title_stands_in_for_a_tag_that_is_not_a_number(tmp_path):
    reader = SdfReader(_tagged(tmp_path / "t.sdf", "REG-0007", "7"))
    assert [r.cid for r in reader] == [7]
    assert reader.skipped_id is None


def test_text_ids_are_named_as_the_reason(tmp_path):
    """ChEMBL, ZINC and vendor files: every record was skipped with no reason given."""
    mol = Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")
    mol.SetProp("_Name", "CHEMBL25")
    path = tmp_path / "chembl.sdf"
    with Chem.SDWriter(str(path)) as writer:
        writer.write(mol)
    reader = SdfReader(path)
    assert list(reader) == []
    assert reader.skipped_id == "CHEMBL25"


def test_title_that_is_not_utf8_does_not_stop_the_file(tmp_path):
    """A Windows-1252 name such as "(±)-Ibuprofen" as the title raised
    UnicodeDecodeError and ended the build, although the ID was in the tag."""
    path = _tagged(tmp_path / "latin.sdf", "3672", "PLACEHOLDER")
    path.write_bytes(
        path.read_bytes().replace(b"PLACEHOLDER", "(\u00b1)-Ibuprofen".encode("cp1252"))
    )
    reader = SdfReader(path)
    assert [r.cid for r in reader] == [3672]

    # With the ID only in the unreadable title, the record is a counted skip.
    mol = Chem.MolFromSmiles("CCO")
    mol.SetProp("_Name", "PLACEHOLDER")
    bare = tmp_path / "bare.sdf"
    with Chem.SDWriter(str(bare)) as writer:
        writer.write(mol)
    bare.write_bytes(bare.read_bytes().replace(b"PLACEHOLDER", b"\xb1 702"))
    reader = SdfReader(bare)
    assert list(reader) == []
    assert (reader.n_failed, reader.skipped_id) == (1, None)


def _many(path, n=400):
    records = [
        (i, "CCO" if i % 2 else "CCN", f"compound {i}", "new") for i in range(1, n + 1)
    ]
    return write_sdf(path, records)


@pytest.mark.parametrize("damage", ["cut", "not gzip", "corrupt"])
def test_damaged_gzip_is_named_as_the_problem(tmp_path, damage):
    """It was reported as "<Boost.Python.function object at 0x...> returned a result
    with an exception set"."""
    from pubchem_finder.core.exceptions import IndexBuildError

    path = _many(tmp_path / "lib.sdf.gz")
    data = bytearray(path.read_bytes())
    if damage == "cut":
        data = data[: len(data) // 2]
    elif damage == "not gzip":
        data = bytearray(b"<html>404 Not Found</html>")
    else:
        data[len(data) // 2] ^= 0xFF
        data[len(data) // 2 + 1] ^= 0xFF
    path.write_bytes(bytes(data))
    reader = SdfReader(path)
    with pytest.raises(IndexBuildError) as err:
        list(reader)
    message = str(err.value)
    assert str(path) in message
    assert "damaged or incomplete" in message
    assert "Boost" not in message


def test_a_stream_that_is_not_deflate_is_named_as_the_problem(tmp_path):
    """zlib.error: the "corrupt" case above fails its CRC, which is an OSError."""
    from pubchem_finder.core.exceptions import IndexBuildError

    path = tmp_path / "lib.sdf.gz"
    path.write_bytes(b"\x1f\x8b\x08\x00\x00\x00\x00\x00\x00\x03\x07" + bytes(32))
    with pytest.raises(IndexBuildError, match="damaged or incomplete"):
        list(SdfReader(path))
