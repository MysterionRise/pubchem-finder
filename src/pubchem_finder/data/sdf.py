"""Streaming reader for PubChem Compound SDF files."""

import gzip
import logging
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import BinaryIO, cast

from rdkit import Chem

logger = logging.getLogger(__name__)

# PubChem renamed its SMILES fields in 2025; accept the new and the old tag names.
SMILES_TAGS = (
    "PUBCHEM_SMILES",
    "PUBCHEM_OPENEYE_ISO_SMILES",
    "PUBCHEM_ISOMERIC_SMILES",
)
CID_TAG = "PUBCHEM_COMPOUND_CID"
INCHIKEY_TAG = "PUBCHEM_IUPAC_INCHIKEY"
IUPAC_NAME_TAG = "PUBCHEM_IUPAC_NAME"
FORMULA_TAG = "PUBCHEM_MOLECULAR_FORMULA"
MW_TAG = "PUBCHEM_MOLECULAR_WEIGHT"
# FPSim2 stores and returns compound IDs as unsigned 32-bit numbers: a larger ID
# would come back from a search as another number (and so as another compound).
MAX_ID = 0xFFFFFFFF


@dataclass(frozen=True, slots=True)
class SdfRecord:
    """One parsed compound."""

    cid: int
    mol: Chem.Mol
    smiles: str
    inchikey: str | None
    iupac_name: str | None
    molecular_formula: str | None
    molecular_weight: float | None


class SdfReader:
    """
    Iterate over the compounds in a (optionally gzipped) SDF file.

    Records that RDKit cannot sanitize, or whose ID is not a whole number from 1
    to :data:`MAX_ID`, are skipped and counted in :attr:`n_failed`;
    :attr:`skipped_id` keeps the first ID that was a number but out of range. The
    whole file is never loaded into memory.

    Example:
        >>> reader = SdfReader("Compound_000000001_000500000.sdf.gz")
        >>> for record in reader: ...
        >>> reader.n_ok, reader.n_failed
    """

    def __init__(self, path: Path | str) -> None:
        self.path = Path(path)
        self.n_ok = 0
        self.n_failed = 0
        self.skipped_id: str | None = None

    def _open(self) -> BinaryIO:
        if self.path.suffix == ".gz":
            return cast(BinaryIO, gzip.open(self.path, "rb"))
        return open(self.path, "rb")

    def __iter__(self) -> Iterator[SdfRecord]:
        with self._open() as stream:
            supplier = Chem.ForwardSDMolSupplier(stream, sanitize=True, removeHs=True)
            for mol in supplier:
                record = _to_record(mol)
                if record is None:
                    self.n_failed += 1
                    if self.skipped_id is None and mol is not None:
                        self.skipped_id = _numeric_id(mol)
                    continue
                self.n_ok += 1
                yield record


def _prop(mol: Chem.Mol, tag: str) -> str | None:
    if mol.HasProp(tag):
        value = str(mol.GetProp(tag)).strip()
        return value or None
    return None


def _ids(mol: Chem.Mol) -> Iterator[str]:
    """Where a record's ID can be: PubChem's CID tag, then the title line."""
    for value in (_prop(mol, CID_TAG), _prop(mol, "_Name")):
        if value:
            yield value


def _cid(mol: Chem.Mol) -> int | None:
    for value in _ids(mol):
        # ASCII digits only (str.isdigit also accepts "²" and "٣"), and short enough
        # to be in range before int() sees it.
        if value.isascii() and value.isdigit() and len(value) <= 10:
            cid = int(value)
            if 1 <= cid <= MAX_ID:
                return cid
    return None


def _numeric_id(mol: Chem.Mol) -> str | None:
    """The ID of a record skipped for its ID alone, shortened for a message."""
    for value in _ids(mol):
        if value.isascii() and value.isdigit():
            return value if len(value) <= 24 else value[:21] + "..."
    return None


def _to_record(mol: Chem.Mol | None) -> SdfRecord | None:
    if mol is None:
        return None
    cid = _cid(mol)
    if cid is None:
        return None
    smiles = next((s for s in (_prop(mol, t) for t in SMILES_TAGS) if s), None)
    if smiles is None:
        smiles = str(Chem.MolToSmiles(mol))
    mw_text = _prop(mol, MW_TAG)
    try:
        mw = float(mw_text) if mw_text else None
    except ValueError:
        mw = None
    return SdfRecord(
        cid=cid,
        mol=mol,
        smiles=smiles,
        inchikey=_prop(mol, INCHIKEY_TAG),
        iupac_name=_prop(mol, IUPAC_NAME_TAG),
        molecular_formula=_prop(mol, FORMULA_TAG),
        molecular_weight=mw,
    )
