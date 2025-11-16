"""MACCS keys fingerprint implementation."""

import numpy as np
from rdkit import Chem
from rdkit.Chem import MACCSkeys

from ..core.fingerprint import FingerprintGenerator


class MACCSFingerprint(FingerprintGenerator):
    """
    MACCS structural keys fingerprint.

    MACCS (Molecular ACCess System) keys are 166 predefined structural keys
    representing the presence or absence of specific substructures.

    The fingerprint is actually 167 bits (0-166), but bit 0 is not used,
    so effectively 166 bits.

    MACCS keys are good for:
    - Diversity analysis
    - Scaffold hopping
    - Fast screening

    Example:
        >>> from rdkit import Chem
        >>> mol = Chem.MolFromSmiles("c1ccccc1")  # benzene
        >>> gen = MACCSFingerprint()
        >>> fp = gen.generate(mol)
        >>> fp.shape
        (167,)
    """

    @property
    def name(self) -> str:
        return "maccs"

    @property
    def size(self) -> int:
        return 167  # MACCS keys are 167 bits (0-166, but 0 is not used)

    def generate(self, mol: Chem.Mol) -> np.ndarray:
        """
        Generate MACCS fingerprint.

        Args:
            mol: RDKit Mol object

        Returns:
            Binary numpy array of length 167
        """
        fp = MACCSkeys.GenMACCSKeys(mol)
        arr = np.zeros((167,), dtype=np.uint8)
        for i in range(167):
            arr[i] = fp[i]
        return arr
