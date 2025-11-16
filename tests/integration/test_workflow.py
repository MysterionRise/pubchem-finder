"""Integration test for complete molecule processing workflow."""

from pubchem_finder.chemistry import MoleculeStandardizer, calculate_properties
from pubchem_finder.core import Molecule
from pubchem_finder.fingerprints import MorganFingerprint, TanimotoSimilarity
from rdkit import Chem


def test_end_to_end_workflow():
    """Test complete workflow from SMILES to similarity comparison."""

    # Test data: aspirin and similar compounds
    molecules_data = [
        ("aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
        ("salicylic_acid", "O=C(O)c1ccccc1O"),
        ("benzoic_acid", "O=C(O)c1ccccc1"),
        ("caffeine", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"),
    ]

    standardizer = MoleculeStandardizer()
    fp_generator = MorganFingerprint(radius=2, n_bits=2048)
    molecules = []

    # Process all molecules
    for name, smiles in molecules_data:
        # Parse
        mol = Chem.MolFromSmiles(smiles)
        assert mol is not None, f"Failed to parse {name}"

        # Standardize
        std_mol = standardizer.standardize(mol)
        assert std_mol is not None, f"Failed to standardize {name}"

        # Calculate properties
        props = calculate_properties(std_mol)

        # Create Molecule object
        molecule = Molecule(
            id=name,
            smiles=smiles,
            canonical_smiles=Chem.MolToSmiles(std_mol),
            **props,
        )

        # Generate fingerprint
        fp = fp_generator.generate(std_mol)
        molecule.add_fingerprint("morgan_2", fp)

        molecules.append(molecule)

    # Verify all molecules processed
    assert len(molecules) == 4

    # Test similarity calculation
    metric = TanimotoSimilarity()
    query_fp = molecules[0].get_fingerprint("morgan_2")  # Aspirin

    for mol in molecules[1:]:
        target_fp = mol.get_fingerprint("morgan_2")
        similarity = metric.calculate(query_fp, target_fp)

        print(f"Aspirin vs {mol.id}: {similarity:.3f}")
        assert 0 <= similarity <= 1

    # Aspirin should be more similar to salicylic acid than to caffeine
    sal_fp = molecules[1].get_fingerprint("morgan_2")
    caf_fp = molecules[3].get_fingerprint("morgan_2")

    sim_salicylic = metric.calculate(query_fp, sal_fp)
    sim_caffeine = metric.calculate(query_fp, caf_fp)

    assert (
        sim_salicylic > sim_caffeine
    ), "Aspirin should be more similar to salicylic acid than caffeine"
