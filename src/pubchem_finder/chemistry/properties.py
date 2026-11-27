"""Molecular descriptors and drug-likeness checks."""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


def unlabelled(mol: Chem.Mol) -> Chem.Mol:
    """Copy of ``mol`` with isotope labels removed, so [2H] counts as plain H.

    Labels such as deuterium are not chemical difficulty or flexibility; RDKit
    would otherwise count them as extra atoms and unknown fragments.
    """
    plain = Chem.Mol(mol)
    for atom in plain.GetAtoms():
        atom.SetIsotope(0)
    return Chem.RemoveHs(plain)


def calculate_properties(mol: Chem.Mol) -> dict[str, float | int]:
    """
    Calculate common molecular properties.

    Keys match the property fields of :class:`pubchem_finder.core.Molecule`,
    so the result can be splatted into its constructor.

    Example:
        >>> mol = Chem.MolFromSmiles("CCO")
        >>> round(calculate_properties(mol)["molecular_weight"], 2)
        46.07
    """
    logp, _molar_refractivity = rdMolDescriptors.CalcCrippenDescriptors(mol)
    return {
        "molecular_weight": float(rdMolDescriptors._CalcMolWt(mol)),
        "logp": float(logp),
        "num_h_donors": int(rdMolDescriptors.CalcNumHBD(mol)),
        "num_h_acceptors": int(rdMolDescriptors.CalcNumHBA(mol)),
        "tpsa": float(rdMolDescriptors.CalcTPSA(mol)),
        "num_rotatable_bonds": int(
            rdMolDescriptors.CalcNumRotatableBonds(unlabelled(mol))
        ),
        "num_aromatic_rings": int(rdMolDescriptors.CalcNumAromaticRings(mol)),
        "num_heavy_atoms": int(mol.GetNumHeavyAtoms()),
    }


def is_druglike(mol: Chem.Mol) -> bool:
    """
    Lipinski's Rule of Five: MW <= 500, LogP <= 5, HBD <= 5, HBA <= 10.
    """
    props = calculate_properties(mol)
    return (
        props["molecular_weight"] <= 500
        and props["logp"] <= 5
        and props["num_h_donors"] <= 5
        and props["num_h_acceptors"] <= 10
    )
