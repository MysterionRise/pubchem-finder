"""Integration tests for molecule processing pipeline."""

import pytest
from pubchem_finder.chemistry import MoleculeStandardizer, calculate_properties
from pubchem_finder.core import Molecule
from pubchem_finder.fingerprints import (
    MACCSFingerprint,
    MorganFingerprint,
    TanimotoSimilarity,
    fingerprint_registry,
)
from rdkit import Chem


class TestMoleculeProcessingPipeline:
    """Test complete molecule processing workflow."""

    def test_full_pipeline_single_molecule(self, aspirin_smiles):
        """Test processing a single molecule through the full pipeline."""
        # 1. Parse SMILES
        mol = Chem.MolFromSmiles(aspirin_smiles)
        assert mol is not None

        # 2. Standardize
        standardizer = MoleculeStandardizer()
        std_mol = standardizer.standardize(mol)
        assert std_mol is not None

        # 3. Calculate properties
        props = calculate_properties(std_mol)
        assert props["molecular_weight"] > 0
        assert "logp" in props

        # 4. Generate canonical SMILES
        canonical_smiles = Chem.MolToSmiles(std_mol)
        assert len(canonical_smiles) > 0

        # 5. Create Molecule object
        molecule = Molecule(
            smiles=aspirin_smiles,
            canonical_smiles=canonical_smiles,
            **props,
        )
        assert molecule.molecular_weight > 0

        # 6. Generate fingerprints
        morgan_gen = MorganFingerprint(radius=2, n_bits=2048)
        morgan_fp = morgan_gen.generate(std_mol)
        molecule.add_fingerprint("morgan_2", morgan_fp)

        maccs_gen = MACCSFingerprint()
        maccs_fp = maccs_gen.generate(std_mol)
        molecule.add_fingerprint("maccs", maccs_fp)

        # Verify everything worked
        assert molecule.has_fingerprint("morgan_2")
        assert molecule.has_fingerprint("maccs")
        assert morgan_fp.shape == (2048,)
        assert maccs_fp.shape == (167,)

    def test_batch_molecule_processing(self):
        """Test processing multiple molecules in batch."""
        smiles_list = ["CCO", "CCCO", "c1ccccc1", "CC(=O)O"]
        molecules = []

        standardizer = MoleculeStandardizer()
        fp_gen = MorganFingerprint(radius=2, n_bits=2048)

        for smiles in smiles_list:
            mol = Chem.MolFromSmiles(smiles)
            std_mol = standardizer.standardize(mol)

            if std_mol:
                props = calculate_properties(std_mol)
                canonical = Chem.MolToSmiles(std_mol)

                molecule = Molecule(
                    smiles=smiles,
                    canonical_smiles=canonical,
                    **props,
                )

                fp = fp_gen.generate(std_mol)
                molecule.add_fingerprint("morgan_2", fp)
                molecules.append(molecule)

        # Verify all molecules processed
        assert len(molecules) == len(smiles_list)
        assert all(m.has_fingerprint("morgan_2") for m in molecules)

    def test_similarity_search_workflow(self):
        """Test similarity search between molecules."""
        # Create query molecule
        query_smiles = "CCO"  # Ethanol
        query_mol = Chem.MolFromSmiles(query_smiles)

        # Create target molecules
        target_smiles = [
            "CCCO",  # Propanol (similar)
            "CC(C)O",  # Isopropanol (similar)
            "c1ccccc1",  # Benzene (dissimilar)
            "C",  # Methane (dissimilar)
        ]

        # Generate fingerprints
        fp_gen = MorganFingerprint(radius=2, n_bits=2048)
        query_fp = fp_gen.generate(query_mol)

        target_fps = [fp_gen.generate(Chem.MolFromSmiles(smi)) for smi in target_smiles]

        # Calculate similarities
        metric = TanimotoSimilarity()
        similarities = [metric.calculate(query_fp, tfp) for tfp in target_fps]

        # Verify results make sense
        assert len(similarities) == len(target_smiles)
        assert all(0 <= sim <= 1 for sim in similarities)

        # Propanol and Isopropanol should be more similar than benzene
        assert similarities[0] > similarities[2]  # CCCO > benzene
        assert similarities[1] > similarities[2]  # CC(C)O > benzene


class TestFingerprintRegistry:
    """Test fingerprint registry integration."""

    def test_registry_with_real_molecules(self, ethanol_mol, benzene_mol):
        """Test using registry to generate multiple fingerprints."""
        # Get all available fingerprints
        available = fingerprint_registry.list_available()
        assert len(available) > 0

        results = {}
        for fp_name in ["morgan_2", "maccs", "ecfp4"]:
            gen = fingerprint_registry.get(fp_name)
            fp = gen.generate(ethanol_mol)
            results[fp_name] = fp

        # Verify all generated
        assert len(results) == 3
        assert results["morgan_2"].shape == (2048,)
        assert results["maccs"].shape == (167,)
        assert results["ecfp4"].shape == (2048,)

    def test_custom_fingerprint_registration(self, ethanol_mol):
        """Test registering and using custom fingerprint."""
        # Create custom Morgan with different parameters
        custom_gen = MorganFingerprint(radius=4, n_bits=4096)
        fingerprint_registry.register("custom_ecfp8", custom_gen)

        # Verify it was registered
        assert "custom_ecfp8" in fingerprint_registry

        # Use it
        gen = fingerprint_registry.get("custom_ecfp8")
        fp = gen.generate(ethanol_mol)

        assert fp.shape == (4096,)


class TestStandardizationEdgeCases:
    """Test standardization with challenging molecules."""

    def test_large_molecule_with_multiple_salts(self):
        """Test standardizing complex molecule with multiple counterions."""
        # Complex molecule with multiple salts
        smiles = "CC(=O)Oc1ccccc1C(=O)O.[Na+].[Na+].[Cl-].[Cl-]"
        standardizer = MoleculeStandardizer()

        result = standardizer.standardize_smiles(smiles)

        # Should keep only the organic molecule
        assert result is not None
        assert "[Na+]" not in result
        assert "[Cl-]" not in result

    def test_charged_molecule_standardization(self):
        """Test standardization handles charges correctly."""
        # Carboxylate anion
        smiles = "CC(=O)[O-]"
        standardizer = MoleculeStandardizer(neutralize=True)

        result = standardizer.standardize_smiles(smiles)
        mol = Chem.MolFromSmiles(result)

        # Should be neutralized
        total_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
        assert total_charge == 0

    def test_aromatic_normalization(self):
        """Test aromatic ring standardization."""
        # Different representations of benzene
        variants = ["c1ccccc1", "C1=CC=CC=C1", "c1cccc1c"]

        standardizer = MoleculeStandardizer()
        results = [standardizer.standardize_smiles(s) for s in variants]

        # All should give same canonical form
        assert len(set(results)) == 1


class TestErrorHandling:
    """Test error handling in integration scenarios."""

    def test_invalid_smiles_in_pipeline(self):
        """Test pipeline handles invalid SMILES gracefully."""
        invalid_smiles = "INVALID_SMILES_XYZ"
        mol = Chem.MolFromSmiles(invalid_smiles)

        # RDKit returns None for invalid SMILES
        assert mol is None

    def test_standardization_failure_handling(self):
        """Test handling of standardization failures."""
        standardizer = MoleculeStandardizer()

        # Test with None input
        result = standardizer.standardize(None)
        assert result is None

        # Test with invalid SMILES
        result = standardizer.standardize_smiles("INVALID")
        assert result is None

    def test_fingerprint_with_invalid_molecule(self):
        """Test fingerprint generation with None molecule."""
        from pubchem_finder.core import FingerprintGenerationError

        gen = MorganFingerprint(radius=2, n_bits=2048)

        # Should raise error for invalid SMILES
        with pytest.raises(FingerprintGenerationError):
            gen.generate_from_smiles("INVALID_SMILES")
