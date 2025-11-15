"""Base fingerprint implementations and similarity metrics."""

import numpy as np

from ..core.fingerprint import SimilarityMetric


class TanimotoSimilarity(SimilarityMetric):
    """
    Tanimoto coefficient for binary fingerprints.

    Also known as Jaccard similarity for binary vectors.
    Formula: |A ∩ B| / |A ∪ B|

    Range: [0, 1] where 1 = identical, 0 = no similarity
    """

    @property
    def name(self) -> str:
        return "tanimoto"

    def calculate(self, fp1: np.ndarray, fp2: np.ndarray) -> float:
        """
        Calculate Tanimoto similarity.

        Args:
            fp1: First fingerprint (binary array)
            fp2: Second fingerprint (binary array)

        Returns:
            Tanimoto coefficient [0, 1]
        """
        intersection = np.sum(fp1 & fp2)
        union = np.sum(fp1 | fp2)
        return float(intersection / union) if union > 0 else 0.0


class DiceSimilarity(SimilarityMetric):
    """
    Dice coefficient for binary fingerprints.

    Formula: 2 * |A ∩ B| / (|A| + |B|)

    Range: [0, 1] where 1 = identical, 0 = no similarity
    """

    @property
    def name(self) -> str:
        return "dice"

    def calculate(self, fp1: np.ndarray, fp2: np.ndarray) -> float:
        """
        Calculate Dice similarity.

        Args:
            fp1: First fingerprint (binary array)
            fp2: Second fingerprint (binary array)

        Returns:
            Dice coefficient [0, 1]
        """
        intersection = np.sum(fp1 & fp2)
        total = np.sum(fp1) + np.sum(fp2)
        return float(2 * intersection / total) if total > 0 else 0.0


class CosineSimilarity(SimilarityMetric):
    """
    Cosine similarity for fingerprints.

    Works with both binary and count fingerprints.
    Formula: (A · B) / (||A|| * ||B||)

    Range: [0, 1] for binary fingerprints
    """

    @property
    def name(self) -> str:
        return "cosine"

    def calculate(self, fp1: np.ndarray, fp2: np.ndarray) -> float:
        """
        Calculate cosine similarity.

        Args:
            fp1: First fingerprint
            fp2: Second fingerprint

        Returns:
            Cosine similarity [0, 1]
        """
        dot_product = np.dot(fp1, fp2)
        norm1 = np.linalg.norm(fp1)
        norm2 = np.linalg.norm(fp2)

        if norm1 == 0 or norm2 == 0:
            return 0.0

        return float(dot_product / (norm1 * norm2))
