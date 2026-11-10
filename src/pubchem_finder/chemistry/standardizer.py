"""Molecule standardization."""

import logging
from dataclasses import dataclass, field

from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

logger = logging.getLogger(__name__)

# Elements that are not metals; every other element is one.
_NON_METALS = frozenset(
    {
        1,
        2,
        5,
        6,
        7,
        8,
        9,
        10,
        14,
        15,
        16,
        17,
        18,
        32,
        33,
        34,
        35,
        36,
        52,
        53,
        54,
        85,
        86,
    }
)
# Groups 1 and 2: counter-ions wherever they are drawn, bonded or not.
_SALT_METALS = frozenset({3, 11, 19, 37, 55, 87, 4, 12, 20, 38, 56, 88})
# Dropped without comment: these are counter-ions far more often than the drug.
_QUIET_METALS = frozenset({11, 19, 12, 20})  # Na, K, Mg, Ca
# RDKit's cleanup draws a sulfoxide as [S+]-[O-]. PubChem, and so the local index,
# draws it S=O, as does everything a user pastes from there.
_CHARGED_SULFOXIDE = Chem.MolFromSmarts("[!O][S+1;X3]([O-1;X1])[!O]")
# Counter-ions that are larger than many of the drugs they are paired with, so
# "keep the largest fragment" would keep them. Neutral forms, as canonical SMILES.
_LARGE_COUNTER_IONS = frozenset(
    str(Chem.MolToSmiles(Chem.MolFromSmiles(smiles)))
    for smiles in (
        "OC(=O)c1cc2ccccc2c(Cc2c(O)c(cc3ccccc23)C(O)=O)c1O",  # pamoic (embonic) acid
    )
)


@dataclass(frozen=True)
class Standardization:
    """A standardized molecule and what standardizing it left out."""

    mol: Chem.Mol
    removed: list[str] = field(default_factory=list)
    note: str | None = None


class MoleculeStandardizer:
    """
    Standardize molecules to a canonical parent form.

    Pipeline (the order matters: the parent fragment is chosen before
    neutralization so counter-ions never influence the charge balance):

    1. ``rdMolStandardize.Cleanup``: remove explicit hydrogens, disconnect
       group 1 and 2 metals, normalize functional groups, reionize, sanitize.
       Sulfoxides stay ``S=O``.
    2. Keep the largest organic fragment, passing over counter-ions known to
       outsize drugs (pamoate). Everything else is reported as removed.
    3. Neutralize charges (optional, default on).

    Two kinds of input have no parent to extract and are returned as given, with
    a note: compounds of a bonded metal other than groups 1 and 2 (cisplatin,
    ferrocene: splitting them at the metal leaves a ligand), and inputs without
    an organic part (sodium chloride).

    Example:
        >>> MoleculeStandardizer().standardize_smiles("CC(=O)[O-].[Na+]")
        'CC(=O)O'
    """

    def __init__(self, neutralize: bool = True) -> None:
        self.neutralize = neutralize
        self._fragment_chooser = rdMolStandardize.LargestFragmentChooser(
            preferOrganic=True
        )
        self._uncharger = rdMolStandardize.Uncharger()

    def standardize(self, mol: Chem.Mol | None) -> Chem.Mol | None:
        """
        Standardize an RDKit molecule.

        Returns:
            The standardized molecule, or ``None`` if the input is ``None`` or
            RDKit cannot process it.
        """
        result = self.standardize_detailed(mol)
        return result.mol if result is not None else None

    def standardize_detailed(self, mol: Chem.Mol | None) -> Standardization | None:
        """Like :meth:`standardize`, also saying which fragments were left out."""
        if mol is None:
            return None
        try:
            return self._standardize(mol)
        except Exception:
            logger.warning("Standardization failed", exc_info=True)
            return None

    def _standardize(self, mol: Chem.Mol) -> Standardization:
        bonded = sorted(
            {
                atom.GetSymbol()
                for atom in mol.GetAtoms()
                if _is_metal(atom)
                and atom.GetAtomicNum() not in _SALT_METALS
                and atom.GetDegree() > 0
            }
        )
        if bonded:
            return Standardization(
                Chem.Mol(mol),
                note=f"Contains bonded {_listed(bonded)}: splitting a metal compound "
                "at the metal would leave a ligand, so it is returned as given.",
            )
        clean = _neutral_sulfoxides(rdMolStandardize.Cleanup(mol))
        fragments = list(Chem.GetMolFrags(clean, asMols=True))
        if not any(_has_carbon(frag) for frag in fragments):
            return Standardization(
                Chem.Mol(mol),
                note="No organic component, so there is no parent to extract: "
                "returned as given.",
            )
        parent = self._fragment_chooser.choose(_without_counter_ions(clean, fragments))
        kept = Chem.MolToSmiles(parent)
        others = [Chem.MolToSmiles(frag) for frag in fragments]
        others.remove(kept)
        if self.neutralize:
            parent = self._uncharger.uncharge(parent)
        parent = _neutral_sulfoxides(parent)
        Chem.SanitizeMol(parent)
        metals = sorted(
            {
                atom.GetSymbol()
                for frag in fragments
                for atom in frag.GetAtoms()
                if _is_metal(atom) and atom.GetAtomicNum() not in _QUIET_METALS
            }
        )
        note = (
            f"Left out {_listed(metals)}: if the metal is the active part, use the "
            "input as it is."
            if metals
            else None
        )
        return Standardization(parent, removed=others, note=note)

    def standardize_smiles(self, smiles: str) -> str | None:
        """
        Standardize a SMILES string and return canonical SMILES.

        Returns:
            Canonical SMILES of the standardized molecule, or ``None`` if the
            SMILES is invalid or standardization fails.
        """
        mol: Chem.Mol | None = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.warning("Invalid SMILES: %r", smiles)
            return None
        std = self.standardize(mol)
        if std is None:
            return None
        return str(Chem.MolToSmiles(std))


def _is_metal(atom: Chem.Atom) -> bool:
    return atom.GetAtomicNum() > 0 and atom.GetAtomicNum() not in _NON_METALS


def _has_carbon(mol: Chem.Mol) -> bool:
    return any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())


def _listed(symbols: list[str]) -> str:
    return ", ".join(symbols)


def _neutral_sulfoxides(mol: Chem.Mol) -> Chem.Mol:
    """``mol`` with every ``[S+]-[O-]`` sulfoxide drawn ``S=O``; stereo is kept."""
    matches = mol.GetSubstructMatches(_CHARGED_SULFOXIDE)
    if not matches:
        return mol
    editable = Chem.RWMol(mol)
    for _, sulfur, oxygen, _ in matches:
        editable.GetAtomWithIdx(sulfur).SetFormalCharge(0)
        editable.GetAtomWithIdx(oxygen).SetFormalCharge(0)
        editable.GetBondBetweenAtoms(sulfur, oxygen).SetBondType(Chem.BondType.DOUBLE)
    result = editable.GetMol()
    Chem.SanitizeMol(result)
    return result


def _without_counter_ions(mol: Chem.Mol, fragments: list[Chem.Mol]) -> Chem.Mol:
    """``mol`` without the fragments that are known large counter-ions, as long as
    an organic fragment is left."""
    uncharger = rdMolStandardize.Uncharger()
    keep = [
        frag
        for frag in fragments
        if Chem.MolToSmiles(uncharger.uncharge(Chem.Mol(frag)))
        not in _LARGE_COUNTER_IONS
    ]
    if len(keep) == len(fragments) or not any(_has_carbon(frag) for frag in keep):
        return mol
    combined = keep[0]
    for frag in keep[1:]:
        combined = Chem.CombineMols(combined, frag)
    return combined
