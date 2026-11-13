"""Streaming reader for PubChem Compound SDF files."""

import gzip
import logging
import zlib
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import BinaryIO, cast

from rdkit import Chem

from ..core.exceptions import IndexBuildError

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
    :attr:`skipped_id` keeps the first ID that could not be used, so the build can
    say why records are missing. The whole file is never loaded into memory.

    Raises:
        IndexBuildError: while iterating, when a gzipped file is cut off or damaged.

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
            # None stands for a record RDKit cannot sanitize.
            supplier: Iterator[Chem.Mol | None] = iter(
                Chem.ForwardSDMolSupplier(stream, sanitize=True, removeHs=True)
            )
            try:
                for mol in supplier:
                    # A record with no atoms has no structure to show or compare.
                    if mol is None or mol.GetNumAtoms() == 0:
                        self.n_failed += 1
                        continue
                    cid, unusable = _cid(mol)
                    if cid is None:
                        self.n_failed += 1
                        if self.skipped_id is None and unusable is not None:
                            self.skipped_id = _shorten(unusable)
                        continue
                    self.n_ok += 1
                    yield _to_record(mol, cid)
            except SystemError as exc:
                # RDKit reads the stream from C++, where an error raised by gzip
                # comes back as a bare SystemError naming a Boost.Python object.
                cause = exc.__cause__ or exc.__context__
                if not isinstance(cause, (EOFError, OSError, zlib.error)):
                    raise
                raise IndexBuildError(
                    f"{self.path} is damaged or incomplete ({cause}); it was read "
                    f"up to record {self.n_ok + self.n_failed:,}"
                ) from exc


def _prop(mol: Chem.Mol, tag: str) -> str | None:
    if mol.HasProp(tag):
        try:
            value = str(mol.GetProp(tag)).strip()
        except UnicodeDecodeError:  # not UTF-8: unreadable, like a missing field
            return None
        return value or None
    return None


def _ids(mol: Chem.Mol) -> Iterator[str]:
    """Where a record's ID can be: PubChem's CID tag, then the title line."""
    for tag in (CID_TAG, "_Name"):
        value = _prop(mol, tag)
        if value:
            yield value


def parse_id(value: str) -> int | None:
    """``value`` as an ID the index can hold (1 to :data:`MAX_ID`), or None."""
    # ASCII digits only: str.isdigit also accepts "²" and "٣".
    if not (value.isascii() and value.isdigit()):
        return None
    # Leading zeros are padding (000000002244); what is left has to be short enough
    # to be in range before int() sees it.
    digits = value.lstrip("0")
    if not 0 < len(digits) <= 10:
        return None
    cid = int(digits)
    return cid if cid <= MAX_ID else None


def _cid(mol: Chem.Mol) -> tuple[int | None, str | None]:
    """A record's ID, or else the first value that could not be used as one.

    The title line stands in when the CID tag is missing or is not a number. A tag
    that is a number the index cannot hold is still the record's ID: the record is
    skipped, not filed under whatever number its title happens to be.
    """
    unusable: str | None = None
    for value in _ids(mol):
        cid = parse_id(value)
        if cid is not None:
            return cid, None
        unusable = unusable or value
        if value.isascii() and value.isdigit():
            break
    return None, unusable


def _shorten(value: str) -> str:
    return value if len(value) <= 24 else value[:21] + "..."


def _to_record(mol: Chem.Mol, cid: int) -> SdfRecord:
    # A tag that does not parse is not SMILES: the index would store and show it,
    # and the ID check could not compare it, so RDKit's own SMILES stands in.
    smiles = next(
        (
            s
            for s in (_prop(mol, t) for t in SMILES_TAGS)
            if s and Chem.MolFromSmiles(s) is not None
        ),
        None,
    )
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
