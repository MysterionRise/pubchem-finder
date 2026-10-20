"""Similarity metrics for binary and count fingerprints."""

from collections.abc import Sequence

import numpy as np
from numpy.typing import NDArray
from rdkit import DataStructs
from rdkit.DataStructs import ExplicitBitVect

from ..core.fingerprint import SimilarityMetric


def _as_bool(fp: NDArray[np.generic]) -> NDArray[np.bool_]:
    return np.asarray(fp).astype(bool)


class TanimotoSimilarity(SimilarityMetric):
    """
    Tanimoto (Jaccard) coefficient for binary fingerprints: |A and B| / |A or B|.

    Non-binary input is treated as binary (any non-zero value is a set bit).
    Range [0, 1]; two empty fingerprints score 0.
    """

    @property
    def name(self) -> str:
        return "tanimoto"

    def calculate(self, fp1: NDArray[np.generic], fp2: NDArray[np.generic]) -> float:
        a, b = self._validate_pair(fp1, fp2)
        a_bits, b_bits = _as_bool(a), _as_bool(b)
        union = int(np.count_nonzero(a_bits | b_bits))
        if union == 0:
            return 0.0
        return int(np.count_nonzero(a_bits & b_bits)) / union

    def batch_calculate(
        self, query_fp: NDArray[np.generic], target_fps: NDArray[np.generic]
    ) -> NDArray[np.float64]:
        q = _as_bool(query_fp)
        targets = _as_bool(target_fps)
        if q.ndim != 1 or targets.ndim != 2 or targets.shape[1] != q.shape[0]:
            raise ValueError(
                "query_fp must be 1-D and target_fps 2-D with matching length"
            )
        inter = np.count_nonzero(targets & q, axis=1)
        union = np.count_nonzero(targets | q, axis=1)
        scores: NDArray[np.float64] = np.divide(
            inter, union, out=np.zeros(len(targets)), where=union > 0
        )
        return scores


class DiceSimilarity(SimilarityMetric):
    """
    Dice coefficient for binary fingerprints: 2|A and B| / (|A| + |B|).

    Range [0, 1]; two empty fingerprints score 0.
    """

    @property
    def name(self) -> str:
        return "dice"

    def calculate(self, fp1: NDArray[np.generic], fp2: NDArray[np.generic]) -> float:
        a, b = self._validate_pair(fp1, fp2)
        a_bits, b_bits = _as_bool(a), _as_bool(b)
        total = int(np.count_nonzero(a_bits)) + int(np.count_nonzero(b_bits))
        if total == 0:
            return 0.0
        return 2 * int(np.count_nonzero(a_bits & b_bits)) / total

    def batch_calculate(
        self, query_fp: NDArray[np.generic], target_fps: NDArray[np.generic]
    ) -> NDArray[np.float64]:
        q = _as_bool(query_fp)
        targets = _as_bool(target_fps)
        if q.ndim != 1 or targets.ndim != 2 or targets.shape[1] != q.shape[0]:
            raise ValueError(
                "query_fp must be 1-D and target_fps 2-D with matching length"
            )
        inter = 2 * np.count_nonzero(targets & q, axis=1)
        total = np.count_nonzero(targets, axis=1) + np.count_nonzero(q)
        scores: NDArray[np.float64] = np.divide(
            inter, total, out=np.zeros(len(targets)), where=total > 0
        )
        return scores


class CosineSimilarity(SimilarityMetric):
    """
    Cosine similarity: (A . B) / (|A| |B|). Works for binary and count vectors.

    Range [0, 1] for non-negative vectors; a zero vector scores 0.
    """

    @property
    def name(self) -> str:
        return "cosine"

    def calculate(self, fp1: NDArray[np.generic], fp2: NDArray[np.generic]) -> float:
        a, b = self._validate_pair(fp1, fp2)
        af = a.astype(float)
        bf = b.astype(float)
        denom = float(np.linalg.norm(af) * np.linalg.norm(bf))
        if denom == 0.0:
            return 0.0
        return float(np.dot(af, bf) / denom)


def bulk_tanimoto(
    query: ExplicitBitVect, targets: Sequence[ExplicitBitVect]
) -> list[float]:
    """Tanimoto similarity of ``query`` to each of ``targets`` using RDKit's C++ loop."""
    return list(DataStructs.BulkTanimotoSimilarity(query, list(targets)))
