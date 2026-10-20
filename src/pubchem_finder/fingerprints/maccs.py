"""MACCS structural keys."""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.DataStructs import ExplicitBitVect

from ..core.fingerprint import FingerprintGenerator

MACCS_SIZE = 167  # RDKit numbers the 166 keys from 1; bit 0 is always unset.


class MACCSFingerprint(FingerprintGenerator):
    """
    MACCS keys: 166 predefined substructure keys (RDKit emits 167 bits).

    Useful for fast screening, diversity analysis and scaffold hopping.

    Example:
        >>> MACCSFingerprint().generate_from_smiles("c1ccccc1").shape
        (167,)
    """

    @property
    def name(self) -> str:
        return "maccs"

    @property
    def size(self) -> int:
        return MACCS_SIZE

    def generate_bitvect(self, mol: Chem.Mol) -> ExplicitBitVect:
        bitvect: ExplicitBitVect = rdMolDescriptors.GetMACCSKeysFingerprint(mol)
        return bitvect
