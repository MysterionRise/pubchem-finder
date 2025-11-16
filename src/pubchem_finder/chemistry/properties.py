"""Calculate molecular descriptors and properties."""

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski


def calculate_properties(mol: Chem.Mol) -> dict[str, float]:
    """
    Calculate common molecular properties.

    These properties can be used for filtering search results
    and applying Lipinski's Rule of Five.

    Args:
        mol: RDKit Mol object

    Returns:
        Dictionary of molecular properties

    Example:
        >>> mol = Chem.MolFromSmiles("CCO")
        >>> props = calculate_properties(mol)
        >>> props['molecular_weight']
        46.07
    """
    return {
        "molecular_weight": Descriptors.MolWt(mol),
        "logp": Crippen.MolLogP(mol),
        "num_h_donors": Lipinski.NumHDonors(mol),
        "num_h_acceptors": Lipinski.NumHAcceptors(mol),
        "tpsa": Descriptors.TPSA(mol),
        "num_rotatable_bonds": Lipinski.NumRotatableBonds(mol),
        "num_aromatic_rings": Lipinski.NumAromaticRings(mol),
        "num_heavy_atoms": Lipinski.HeavyAtomCount(mol),
    }


def is_druglike(mol: Chem.Mol) -> bool:
    """
    Check if molecule satisfies Lipinski's Rule of Five.

    Rules:
    - Molecular weight <= 500 Da
    - LogP <= 5
    - Number of H-bond donors <= 5
    - Number of H-bond acceptors <= 10

    Args:
        mol: RDKit Mol object

    Returns:
        True if molecule is drug-like
    """
    props = calculate_properties(mol)
    return (
        props["molecular_weight"] <= 500
        and props["logp"] <= 5
        and props["num_h_donors"] <= 5
        and props["num_h_acceptors"] <= 10
    )
