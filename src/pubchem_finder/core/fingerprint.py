"""Abstract base classes for fingerprint generators and similarity metrics."""

from abc import ABC, abstractmethod

import numpy as np
from numpy.typing import NDArray
from rdkit import Chem, DataStructs
from rdkit.DataStructs import ExplicitBitVect

from .exceptions import FingerprintGenerationError

BitArray = NDArray[np.uint8]


class FingerprintGenerator(ABC):
    """
    Base class for all fingerprint generators.

    Subclasses implement :meth:`generate_bitvect`, which returns RDKit's native
    ``ExplicitBitVect``. :meth:`generate` converts it to a dense ``uint8`` numpy
    array for callers that want vector maths.
    """

    @property
    @abstractmethod
    def name(self) -> str:
        """Unique identifier for this fingerprint configuration."""

    @property
    @abstractmethod
    def size(self) -> int:
        """Number of bits in the fingerprint."""

    @abstractmethod
    def generate_bitvect(self, mol: Chem.Mol) -> ExplicitBitVect:
        """Generate the fingerprint as an RDKit ``ExplicitBitVect``."""

    def generate(self, mol: Chem.Mol) -> BitArray:
        """
        Generate the fingerprint as a dense ``uint8`` array of length ``size``.

        Raises:
            FingerprintGenerationError: if ``mol`` is ``None`` or RDKit fails.
        """
        if mol is None:
            raise FingerprintGenerationError("Cannot generate a fingerprint from None")
        try:
            bitvect = self.generate_bitvect(mol)
        except Exception as exc:
            raise FingerprintGenerationError(
                f"{self.name}: fingerprint generation failed: {exc}"
            ) from exc
        arr: BitArray = np.zeros((self.size,), dtype=np.uint8)
        DataStructs.ConvertToNumpyArray(bitvect, arr)
        return arr

    def generate_bitvect_from_smiles(self, smiles: str) -> ExplicitBitVect:
        """Parse ``smiles`` and return its fingerprint as an ``ExplicitBitVect``."""
        return self.generate_bitvect(self._mol_from_smiles(smiles))

    def generate_from_smiles(self, smiles: str) -> BitArray:
        """Parse ``smiles`` and return its fingerprint as a numpy array."""
        return self.generate(self._mol_from_smiles(smiles))

    @staticmethod
    def _mol_from_smiles(smiles: str) -> Chem.Mol:
        mol: Chem.Mol | None = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise FingerprintGenerationError(f"Invalid SMILES: {smiles!r}")
        return mol

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(name={self.name}, size={self.size})"


class SimilarityMetric(ABC):
    """Base class for similarity metrics over numpy fingerprints."""

    @property
    @abstractmethod
    def name(self) -> str:
        """Metric identifier."""

    @abstractmethod
    def calculate(self, fp1: NDArray[np.generic], fp2: NDArray[np.generic]) -> float:
        """Similarity between two 1-D fingerprints of equal length."""

    def batch_calculate(
        self, query_fp: NDArray[np.generic], target_fps: NDArray[np.generic]
    ) -> NDArray[np.float64]:
        """
        Similarity between one query and many targets.

        Args:
            query_fp: 1-D array of length ``n``.
            target_fps: 2-D array of shape ``(m, n)``.

        Returns:
            1-D array of ``m`` similarity scores.
        """
        targets = np.asarray(target_fps)
        if targets.ndim != 2:
            raise ValueError(f"target_fps must be 2-D, got shape {targets.shape}")
        return np.array([self.calculate(query_fp, fp) for fp in targets], dtype=float)

    @staticmethod
    def _validate_pair(
        fp1: NDArray[np.generic], fp2: NDArray[np.generic]
    ) -> tuple[NDArray[np.generic], NDArray[np.generic]]:
        a = np.asarray(fp1)
        b = np.asarray(fp2)
        if a.ndim != 1 or b.ndim != 1:
            raise ValueError("Fingerprints must be 1-D arrays")
        if a.shape != b.shape:
            raise ValueError(
                f"Fingerprint lengths differ: {a.shape[0]} vs {b.shape[0]}"
            )
        return a, b

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(name={self.name})"
