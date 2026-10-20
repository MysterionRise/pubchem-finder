"""Tests for fingerprint generators."""

import numpy as np
import pytest

from pubchem_finder.core import FingerprintGenerationError
from pubchem_finder.fingerprints import (
    DiceSimilarity,
    MACCSFingerprint,
    MorganFingerprint,
    TanimotoSimilarity,
    fingerprint_registry,
)


class TestMorganFingerprint:
    """Tests for Morgan fingerprint generator."""

    def test_generate_correct_shape(self, ethanol_mol):
        """Test Morgan fingerprint generates correct shape."""
        generator = MorganFingerprint(radius=2, n_bits=2048)
        fp = generator.generate(ethanol_mol)

        assert fp.shape == (2048,)
        assert fp.dtype == np.uint8

    def test_generate_deterministic(self, ethanol_mol):
        """Test fingerprint is deterministic."""
        generator = MorganFingerprint(radius=2, n_bits=2048)

        fp1 = generator.generate(ethanol_mol)
        fp2 = generator.generate(ethanol_mol)

        np.testing.assert_array_equal(fp1, fp2)

    def test_generate_from_smiles(self, ethanol_smiles):
        """Test generating from SMILES string."""
        generator = MorganFingerprint(radius=2, n_bits=2048)
        fp = generator.generate_from_smiles(ethanol_smiles)

        assert fp.shape == (2048,)

    def test_generate_from_invalid_smiles(self, invalid_smiles):
        """Test error handling for invalid SMILES."""
        generator = MorganFingerprint(radius=2, n_bits=2048)

        with pytest.raises(FingerprintGenerationError):
            generator.generate_from_smiles(invalid_smiles)

    def test_different_radius(self, aspirin_mol):
        """Test fingerprints with different radius."""
        gen_r2 = MorganFingerprint(radius=2, n_bits=2048)
        gen_r3 = MorganFingerprint(radius=3, n_bits=2048)

        fp_r2 = gen_r2.generate(aspirin_mol)
        fp_r3 = gen_r3.generate(aspirin_mol)

        # Different radius should give different fingerprints for larger molecules
        assert not np.array_equal(fp_r2, fp_r3)

    def test_feature_based_morgan(self, ethanol_mol):
        """Test feature-based Morgan (FCFP)."""
        gen_ecfp = MorganFingerprint(radius=2, use_features=False)
        gen_fcfp = MorganFingerprint(radius=2, use_features=True)

        fp_ecfp = gen_ecfp.generate(ethanol_mol)
        fp_fcfp = gen_fcfp.generate(ethanol_mol)

        # ECFP and FCFP should be different
        assert not np.array_equal(fp_ecfp, fp_fcfp)

    def test_name_property(self):
        """Test name property."""
        gen = MorganFingerprint(radius=2, n_bits=2048)
        assert "morgan" in gen.name
        assert "r2" in gen.name
        assert "2048" in gen.name

    def test_size_property(self):
        """Test size property."""
        gen = MorganFingerprint(radius=2, n_bits=2048)
        assert gen.size == 2048


class TestMACCSFingerprint:
    """Tests for MACCS keys fingerprint."""

    def test_generate_correct_shape(self, ethanol_mol):
        """Test MACCS fingerprint is 167 bits."""
        generator = MACCSFingerprint()
        fp = generator.generate(ethanol_mol)

        assert fp.shape == (167,)
        assert fp.dtype == np.uint8

    def test_generate_deterministic(self, benzene_mol):
        """Test MACCS fingerprint is deterministic."""
        generator = MACCSFingerprint()

        fp1 = generator.generate(benzene_mol)
        fp2 = generator.generate(benzene_mol)

        np.testing.assert_array_equal(fp1, fp2)

    def test_generate_from_smiles(self, benzene_smiles):
        """Test generating MACCS from SMILES."""
        generator = MACCSFingerprint()
        fp = generator.generate_from_smiles(benzene_smiles)

        assert fp.shape == (167,)

    def test_name_property(self):
        """Test MACCS name."""
        gen = MACCSFingerprint()
        assert gen.name == "maccs"

    def test_size_property(self):
        """Test MACCS size."""
        gen = MACCSFingerprint()
        assert gen.size == 167


class TestSimilarityMetrics:
    """Tests for similarity metrics."""

    def test_tanimoto_identical(self):
        """Test Tanimoto similarity of identical fingerprints."""
        fp1 = np.array([1, 0, 1, 1, 0])
        fp2 = np.array([1, 0, 1, 1, 0])

        metric = TanimotoSimilarity()
        sim = metric.calculate(fp1, fp2)

        assert sim == 1.0

    def test_tanimoto_no_similarity(self):
        """Test Tanimoto similarity with no overlap."""
        fp1 = np.array([1, 0, 0, 0, 0])
        fp2 = np.array([0, 1, 1, 1, 1])

        metric = TanimotoSimilarity()
        sim = metric.calculate(fp1, fp2)

        assert sim == 0.0

    def test_tanimoto_partial_similarity(self):
        """Test Tanimoto similarity with partial overlap."""
        fp1 = np.array([1, 0, 1, 1, 0])
        fp2 = np.array([1, 1, 1, 0, 0])

        metric = TanimotoSimilarity()
        sim = metric.calculate(fp1, fp2)

        # Intersection: 2 bits, Union: 4 bits -> 2/4 = 0.5
        assert sim == 0.5

    def test_dice_identical(self):
        """Test Dice similarity of identical fingerprints."""
        fp1 = np.array([1, 0, 1, 1, 0])
        fp2 = np.array([1, 0, 1, 1, 0])

        metric = DiceSimilarity()
        sim = metric.calculate(fp1, fp2)

        assert sim == 1.0

    def test_dice_no_similarity(self):
        """Test Dice similarity with no overlap."""
        fp1 = np.array([1, 0, 0, 0, 0])
        fp2 = np.array([0, 1, 1, 1, 1])

        metric = DiceSimilarity()
        sim = metric.calculate(fp1, fp2)

        assert sim == 0.0

    def test_tanimoto_name(self):
        """Test Tanimoto metric name."""
        metric = TanimotoSimilarity()
        assert metric.name == "tanimoto"

    def test_dice_name(self):
        """Test Dice metric name."""
        metric = DiceSimilarity()
        assert metric.name == "dice"


class TestFingerprintRegistry:
    """Tests for fingerprint registry."""

    def test_get_morgan(self):
        """Test getting Morgan fingerprint from registry."""
        gen = fingerprint_registry.get("morgan_2")
        assert isinstance(gen, MorganFingerprint)
        assert gen.radius == 2

    def test_get_maccs(self):
        """Test getting MACCS from registry."""
        gen = fingerprint_registry.get("maccs")
        assert isinstance(gen, MACCSFingerprint)

    def test_get_ecfp4_alias(self):
        """Test ECFP4 alias."""
        gen = fingerprint_registry.get("ecfp4")
        assert isinstance(gen, MorganFingerprint)
        assert gen.radius == 2

    def test_get_nonexistent(self):
        """Test error for non-existent fingerprint."""
        with pytest.raises(KeyError):
            fingerprint_registry.get("nonexistent_fp")

    def test_list_available(self):
        """Test listing available fingerprints."""
        available = fingerprint_registry.list_available()
        assert "morgan_2" in available
        assert "maccs" in available
        assert "ecfp4" in available

    def test_contains(self):
        """Test __contains__ method."""
        assert "morgan_2" in fingerprint_registry
        assert "maccs" in fingerprint_registry
        assert "nonexistent" not in fingerprint_registry

    def test_register_custom(self):
        """Test registering custom fingerprint."""
        custom_gen = MorganFingerprint(radius=4, n_bits=1024)
        fingerprint_registry.register("custom_morgan", custom_gen)

        gen = fingerprint_registry.get("custom_morgan")
        assert gen.radius == 4
        assert gen.size == 1024
