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

    Records that RDKit cannot sanitize, or that carry no integer CID, are
    skipped and counted in :attr:`n_failed`. The whole file is never loaded
    into memory.

    Example:
        >>> reader = SdfReader("Compound_000000001_000500000.sdf.gz")
        >>> for record in reader: ...
        >>> reader.n_ok, reader.n_failed
    """

    def __init__(self, path: Path | str) -> None:
        self.path = Path(path)
        self.n_ok = 0
        self.n_failed = 0

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
                    continue
                self.n_ok += 1
                yield record


def _prop(mol: Chem.Mol, tag: str) -> str | None:
    if mol.HasProp(tag):
        value = str(mol.GetProp(tag)).strip()
        return value or None
    return None


def _cid(mol: Chem.Mol) -> int | None:
    for value in (_prop(mol, CID_TAG), _prop(mol, "_Name")):
        if value and value.isdigit():
            return int(value)
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
