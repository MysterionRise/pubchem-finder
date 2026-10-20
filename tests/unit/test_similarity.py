"""Tests for similarity metrics and bulk helpers."""

import numpy as np
import pytest
from rdkit.Chem import rdFingerprintGenerator

from pubchem_finder.core import FingerprintGenerationError
from pubchem_finder.fingerprints import (
    CosineSimilarity,
    DiceSimilarity,
    MorganFingerprint,
    TanimotoSimilarity,
    bulk_tanimoto,
)


@pytest.mark.parametrize("metric", [TanimotoSimilarity(), DiceSimilarity()])
def test_binary_metrics_accept_float_and_int_arrays(metric):
    a = np.array([1.0, 0.0, 1.0, 1.0])
    b = np.array([1, 1, 1, 0], dtype=np.int64)
    assert metric.calculate(a, b) == metric.calculate(a.astype(np.uint8), b)


@pytest.mark.parametrize(
    "metric", [TanimotoSimilarity(), DiceSimilarity(), CosineSimilarity()]
)
def test_metrics_reject_mismatched_lengths(metric):
    with pytest.raises(ValueError, match="lengths differ"):
        metric.calculate(np.array([1, 0, 1]), np.array([1, 0]))


def test_cosine_identical_and_orthogonal():
    metric = CosineSimilarity()
    assert metric.name == "cosine"
    assert metric.calculate(np.array([1, 0, 1]), np.array([1, 0, 1])) == pytest.approx(
        1
    )
    assert metric.calculate(np.array([1, 0]), np.array([0, 1])) == 0.0
    assert metric.calculate(np.array([0, 0]), np.array([1, 1])) == 0.0


def test_cosine_counts():
    metric = CosineSimilarity()
    assert metric.calculate(np.array([2, 0]), np.array([4, 0])) == pytest.approx(1.0)


@pytest.mark.parametrize("metric", [TanimotoSimilarity(), DiceSimilarity()])
def test_batch_matches_pairwise(metric):
    rng = np.random.default_rng(0)
    query = rng.integers(0, 2, size=64, dtype=np.uint8)
    targets = rng.integers(0, 2, size=(10, 64), dtype=np.uint8)
    targets[3] = 0  # an empty fingerprint must not divide by zero
    batch = metric.batch_calculate(query, targets)
    assert batch.shape == (10,)
    expected = [metric.calculate(query, t) for t in targets]
    np.testing.assert_allclose(batch, expected)


def test_batch_rejects_bad_shapes():
    with pytest.raises(ValueError):
        TanimotoSimilarity().batch_calculate(np.array([1, 0]), np.array([1, 0]))
    with pytest.raises(ValueError):
        CosineSimilarity().batch_calculate(np.array([1, 0]), np.array([1, 0]))


def test_bulk_tanimoto_matches_numpy_metric(aspirin_smiles):
    gen = MorganFingerprint(radius=2, n_bits=1024)
    smiles = [aspirin_smiles, "O=C(O)c1ccccc1O", "CCO", "c1ccccc1"]
    bitvects = [gen.generate_bitvect_from_smiles(s) for s in smiles]
    arrays = np.array([gen.generate_from_smiles(s) for s in smiles])

    rdkit_scores = bulk_tanimoto(bitvects[0], bitvects)
    numpy_scores = TanimotoSimilarity().batch_calculate(arrays[0], arrays)

    assert rdkit_scores[0] == pytest.approx(1.0)
    np.testing.assert_allclose(rdkit_scores, numpy_scores)


def test_generate_none_raises():
    with pytest.raises(FingerprintGenerationError):
        MorganFingerprint().generate(None)


def test_morgan_matches_rdkit_generator(aspirin_mol):
    """Our ECFP4/1024 bits must equal RDKit's generator with the same params.

    FPSim2 builds its index with ``rdFingerprintGenerator`` using
    ``{"radius": 2, "fpSize": 1024}``; query and index fingerprints must agree.
    """
    ours = MorganFingerprint(radius=2, n_bits=1024).generate_bitvect(aspirin_mol)
    reference = rdFingerprintGenerator.GetMorganGenerator(
        radius=2, fpSize=1024
    ).GetFingerprint(aspirin_mol)
    assert ours.ToBitString() == reference.ToBitString()


def test_morgan_rejects_bad_params():
    with pytest.raises(ValueError):
        MorganFingerprint(radius=-1)
    with pytest.raises(ValueError):
        MorganFingerprint(n_bits=0)
