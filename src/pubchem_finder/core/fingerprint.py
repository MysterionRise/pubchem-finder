"""Abstract base classes for fingerprints and similarity metrics."""

from abc import ABC, abstractmethod

import numpy as np
from rdkit import Chem

from .exceptions import FingerprintGenerationError


class FingerprintGenerator(ABC):
    """
    Abstract base class for all fingerprint generators.

    This interface ensures all fingerprint implementations are interchangeable.
    Future fingerprint types (including ML-based) will implement this interface.
    """

    @property
    @abstractmethod
    def name(self) -> str:
        """Unique identifier for this fingerprint type."""
        pass

    @property
    @abstractmethod
    def size(self) -> int:
        """Dimensionality of the fingerprint."""
        pass

    @abstractmethod
    def generate(self, mol: Chem.Mol) -> np.ndarray:
        """
        Generate fingerprint from RDKit molecule.

        Args:
            mol: RDKit Mol object

        Returns:
            numpy array representing the fingerprint

        Raises:
            FingerprintGenerationError: If fingerprint cannot be generated
        """
        pass

    def generate_from_smiles(self, smiles: str) -> np.ndarray:
        """
        Convenience method to generate fingerprint from SMILES.

        Args:
            smiles: SMILES string

        Returns:
            numpy array representing the fingerprint

        Raises:
            FingerprintGenerationError: If SMILES is invalid or fingerprint
                cannot be generated
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise FingerprintGenerationError(f"Invalid SMILES: {smiles}")
        return self.generate(mol)

    def __repr__(self) -> str:
        """String representation."""
        return f"{self.__class__.__name__}(name={self.name}, size={self.size})"


class SimilarityMetric(ABC):
    """
    Abstract base for similarity metrics.

    Extensible design allows adding new metrics without changing search logic.
    """

    @property
    @abstractmethod
    def name(self) -> str:
        """Metric identifier."""
        pass

    @abstractmethod
    def calculate(self, fp1: np.ndarray, fp2: np.ndarray) -> float:
        """
        Calculate similarity between two fingerprints.

        Args:
            fp1: First fingerprint
            fp2: Second fingerprint

        Returns:
            Similarity score (higher = more similar, range depends on metric)
        """
        pass

    def batch_calculate(
        self, query_fp: np.ndarray, target_fps: np.ndarray
    ) -> np.ndarray:
        """
        Calculate similarity between one query and multiple targets.

        Args:
            query_fp: Query fingerprint (1D array)
            target_fps: Target fingerprints (2D array, shape: [n_targets, fp_size])

        Returns:
            Array of similarity scores
        """
        return np.array([self.calculate(query_fp, fp) for fp in target_fps])

    def __repr__(self) -> str:
        """String representation."""
        return f"{self.__class__.__name__}(name={self.name})"
