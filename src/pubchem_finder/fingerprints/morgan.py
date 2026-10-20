"""Morgan (ECFP/FCFP) fingerprints via RDKit's fingerprint generators."""

from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from rdkit.DataStructs import ExplicitBitVect

from ..core.fingerprint import FingerprintGenerator


class MorganFingerprint(FingerprintGenerator):
    """
    Morgan circular fingerprint (ECFP, or FCFP with ``use_features=True``).

    Radius 2 corresponds to ECFP4, radius 3 to ECFP6.

    Args:
        radius: Circular environment radius (default 2, i.e. ECFP4).
        n_bits: Folded fingerprint length (default 2048).
        use_features: Use pharmacophoric feature invariants (FCFP) instead of
            atom-type invariants (ECFP).

    Example:
        >>> gen = MorganFingerprint(radius=2, n_bits=1024)
        >>> gen.generate_from_smiles("CCO").shape
        (1024,)
    """

    def __init__(
        self, radius: int = 2, n_bits: int = 2048, use_features: bool = False
    ) -> None:
        if radius < 0:
            raise ValueError("radius must be >= 0")
        if n_bits <= 0:
            raise ValueError("n_bits must be > 0")
        self.radius = radius
        self.n_bits = n_bits
        self.use_features = use_features
        invariants = (
            rdFingerprintGenerator.GetMorganFeatureAtomInvGen()
            if use_features
            else None
        )
        self._generator = rdFingerprintGenerator.GetMorganGenerator(
            radius=radius, fpSize=n_bits, atomInvariantsGenerator=invariants
        )

    @property
    def name(self) -> str:
        kind = "fcfp" if self.use_features else "ecfp"
        return f"morgan_{kind}_r{self.radius}_{self.n_bits}"

    @property
    def size(self) -> int:
        return self.n_bits

    def generate_bitvect(self, mol: Chem.Mol) -> ExplicitBitVect:
        bitvect: ExplicitBitVect = self._generator.GetFingerprint(mol)
        return bitvect
