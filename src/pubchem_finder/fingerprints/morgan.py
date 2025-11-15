"""Morgan (ECFP) fingerprint implementation."""

import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem

from ..core.fingerprint import FingerprintGenerator


class MorganFingerprint(FingerprintGenerator):
    """
    Morgan fingerprint generator (aka ECFP - Extended Connectivity Fingerprint).

    Morgan fingerprints encode circular atom environments and are widely used
    for molecular similarity searching.

    Radius equivalence:
    - radius=2 corresponds to ECFP4 (diameter=4)
    - radius=3 corresponds to ECFP6 (diameter=6)

    Args:
        radius: Fingerprint radius (default: 2 for ECFP4)
        n_bits: Number of bits in fingerprint (default: 2048)
        use_features: Use feature-based invariants instead of atom types
                      (FCFP vs ECFP, default: False)

    Example:
        >>> from rdkit import Chem
        >>> mol = Chem.MolFromSmiles("CCO")
        >>> gen = MorganFingerprint(radius=2, n_bits=2048)
        >>> fp = gen.generate(mol)
        >>> fp.shape
        (2048,)
    """

    def __init__(
        self, radius: int = 2, n_bits: int = 2048, use_features: bool = False
    ):
        self.radius = radius
        self.n_bits = n_bits
        self.use_features = use_features

    @property
    def name(self) -> str:
        feature_str = "fcfp" if self.use_features else "ecfp"
        return f"morgan_{feature_str}_r{self.radius}_{self.n_bits}"

    @property
    def size(self) -> int:
        return self.n_bits

    def generate(self, mol: Chem.Mol) -> np.ndarray:
        """
        Generate Morgan fingerprint.

        Args:
            mol: RDKit Mol object

        Returns:
            Binary numpy array of length n_bits
        """
        fp = AllChem.GetMorganFingerprintAsBitVect(
            mol, radius=self.radius, nBits=self.n_bits, useFeatures=self.use_features
        )

        # Convert to numpy array
        arr = np.zeros((self.n_bits,), dtype=np.uint8)
        DataStructs.ConvertToNumpyArray(fp, arr)
        return arr
