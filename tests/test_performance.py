"""Performance and benchmark tests."""

import time

import numpy as np
import pytest
from rdkit import Chem

from pubchem_finder.chemistry import MoleculeStandardizer, calculate_properties
from pubchem_finder.fingerprints import (
    MACCSFingerprint,
    MorganFingerprint,
    TanimotoSimilarity,
)


class TestPerformance:
    """Performance benchmarks for core operations."""

    @pytest.fixture
    def molecule_batch(self) -> list[Chem.Mol]:
        """Generate batch of molecules for testing."""
        smiles_list = [
            "CCO",
            "CCCO",
            "c1ccccc1",
            "CC(=O)O",
            "CC(C)O",
            "c1ccc2c(c1)cccc2",
            "CC(=O)Oc1ccccc1C(=O)O",
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        ] * 100  # 800 molecules

        return [Chem.MolFromSmiles(smi) for smi in smiles_list]

    def test_standardization_performance(self, molecule_batch, benchmark):
        """Benchmark molecule standardization."""
        standardizer = MoleculeStandardizer()

        def standardize_batch():
            return [standardizer.standardize(mol) for mol in molecule_batch[:100]]

        if benchmark:
            result = benchmark(standardize_batch)
        else:
            start = time.time()
            result = standardize_batch()
            elapsed = time.time() - start
            throughput = len(result) / elapsed
            print(f"\nStandardization: {throughput:.0f} molecules/sec")

        assert len(result) == 100

    def test_property_calculation_performance(self, molecule_batch, benchmark):
        """Benchmark property calculation."""

        def calculate_batch_properties():
            return [calculate_properties(mol) for mol in molecule_batch[:100]]

        if benchmark:
            result = benchmark(calculate_batch_properties)
        else:
            start = time.time()
            result = calculate_batch_properties()
            elapsed = time.time() - start
            throughput = len(result) / elapsed
            print(f"\nProperty calculation: {throughput:.0f} molecules/sec")

        assert len(result) == 100

    def test_morgan_fingerprint_performance(self, molecule_batch, benchmark):
        """Benchmark Morgan fingerprint generation."""
        generator = MorganFingerprint(radius=2, n_bits=2048)

        def generate_batch_fingerprints():
            return [generator.generate(mol) for mol in molecule_batch[:100]]

        if benchmark:
            result = benchmark(generate_batch_fingerprints)
        else:
            start = time.time()
            result = generate_batch_fingerprints()
            elapsed = time.time() - start
            throughput = len(result) / elapsed
            print(f"\nMorgan fingerprint: {throughput:.0f} molecules/sec")

        assert len(result) == 100

    def test_maccs_fingerprint_performance(self, molecule_batch, benchmark):
        """Benchmark MACCS fingerprint generation."""
        generator = MACCSFingerprint()

        def generate_batch_fingerprints():
            return [generator.generate(mol) for mol in molecule_batch[:100]]

        if benchmark:
            result = benchmark(generate_batch_fingerprints)
        else:
            start = time.time()
            result = generate_batch_fingerprints()
            elapsed = time.time() - start
            throughput = len(result) / elapsed
            print(f"\nMACCS fingerprint: {throughput:.0f} molecules/sec")

        assert len(result) == 100

    def test_tanimoto_similarity_performance(self, benchmark):
        """Benchmark Tanimoto similarity calculation."""
        # Generate random fingerprints
        np.random.seed(42)
        query_fp = np.random.randint(0, 2, size=2048, dtype=np.uint8)
        target_fps = np.random.randint(0, 2, size=(1000, 2048), dtype=np.uint8)

        metric = TanimotoSimilarity()

        def calculate_batch_similarities():
            return [metric.calculate(query_fp, tfp) for tfp in target_fps]

        if benchmark:
            result = benchmark(calculate_batch_similarities)
        else:
            start = time.time()
            result = calculate_batch_similarities()
            elapsed = time.time() - start
            throughput = len(result) / elapsed
            print(f"\nTanimoto similarity: {throughput:.0f} comparisons/sec")

        assert len(result) == 1000


@pytest.fixture
def benchmark():
    """Optional pytest-benchmark fixture."""
    import importlib.util

    if importlib.util.find_spec("pytest_benchmark") is not None:
        return pytest.mark.benchmark
    return None


class TestScalability:
    """Test behavior with large datasets."""

    def test_large_batch_processing(self):
        """Test processing 1000 molecules."""
        # Generate simple molecules - valid SMILES strings
        base_smiles = ["C", "CC", "CCC", "CCCC", "c1ccccc1", "CCO", "CCCO", "CC(C)C", "CCCCC", "CCCCCC"]
        smiles_list = base_smiles * 100  # 1000 molecules

        standardizer = MoleculeStandardizer()
        fp_gen = MorganFingerprint(radius=2, n_bits=2048)

        processed = 0
        for smiles in smiles_list:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                std_mol = standardizer.standardize(mol)
                if std_mol:
                    fp = fp_gen.generate(std_mol)
                    assert fp.shape == (2048,)
                    processed += 1

        assert processed == len(smiles_list)

    def test_memory_efficiency(self):
        """Test that fingerprints don't cause memory issues."""
        import gc

        fp_gen = MorganFingerprint(radius=2, n_bits=2048)
        fingerprints = []

        # Generate 1000 fingerprints
        for i in range(1000):
            mol = Chem.MolFromSmiles("C" * (i % 10 + 1))
            fp = fp_gen.generate(mol)
            fingerprints.append(fp)

            # Periodic garbage collection
            if i % 100 == 0:
                gc.collect()

        assert len(fingerprints) == 1000
        # Each fingerprint should be small
        assert fingerprints[0].nbytes < 5000  # < 5KB


class TestAccuracy:
    """Test accuracy of fingerprints and similarity calculations."""

    def test_identical_molecules_similarity(self):
        """Test that identical molecules have similarity = 1.0."""
        mol = Chem.MolFromSmiles("CCO")
        fp_gen = MorganFingerprint(radius=2, n_bits=2048)

        fp1 = fp_gen.generate(mol)
        fp2 = fp_gen.generate(mol)

        metric = TanimotoSimilarity()
        similarity = metric.calculate(fp1, fp2)

        assert similarity == 1.0

    def test_known_similarity_values(self):
        """Test similarity calculations with known reference values."""
        # Ethanol vs Methanol should have moderate similarity
        ethanol = Chem.MolFromSmiles("CCO")
        methanol = Chem.MolFromSmiles("CO")

        fp_gen = MorganFingerprint(radius=2, n_bits=2048)
        fp_ethanol = fp_gen.generate(ethanol)
        fp_methanol = fp_gen.generate(methanol)

        metric = TanimotoSimilarity()
        similarity = metric.calculate(fp_ethanol, fp_methanol)

        # Should have reasonable similarity (small alcohols are similar)
        # Tanimoto similarity for ethanol/methanol is ~0.28-0.30
        assert 0.2 < similarity < 1.0
        assert similarity > 0

    def test_dissimilar_molecules(self):
        """Test that very different molecules have low similarity."""
        # Small aliphatic vs large aromatic
        methane = Chem.MolFromSmiles("C")
        naphthalene = Chem.MolFromSmiles("c1ccc2c(c1)cccc2")

        fp_gen = MorganFingerprint(radius=2, n_bits=2048)
        fp1 = fp_gen.generate(methane)
        fp2 = fp_gen.generate(naphthalene)

        metric = TanimotoSimilarity()
        similarity = metric.calculate(fp1, fp2)

        # Should have low similarity
        assert similarity < 0.3
