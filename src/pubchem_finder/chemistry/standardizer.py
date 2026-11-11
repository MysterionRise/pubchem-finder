"""Molecule standardization."""

import logging
from dataclasses import dataclass, field

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
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
# Groups 1 and 2: counter-ions wherever they are drawn, unless bonded to carbon.
_SALT_METALS = frozenset({3, 11, 19, 37, 55, 87, 4, 12, 20, 38, 56, 88})
# Dropped without comment: these are counter-ions far more often than the drug.
_QUIET_METALS = frozenset({11, 19, 12, 20})  # Na, K, Mg, Ca
# RDKit's cleanup draws a sulfoxide as [S+]-[O-]. PubChem, and so the local index,
# draws it S=O, as does everything a user pastes from there.
_CHARGED_SULFOXIDE = Chem.MolFromSmarts("[!O][S+1;X3]([O-1;X1])[!O]")


def _flat(mol: Chem.Mol) -> str:
    """A fragment as its neutral form without stereo: how partners are compared."""
    neutral = rdMolStandardize.Uncharger().uncharge(Chem.Mol(mol))
    return str(Chem.MolToSmiles(neutral, isomericSmiles=False))


# Salt formers and solvents: what a drug is formulated or crystallized with. Any of
# them can outsize the drug it is paired with (pamoic acid has 28 heavy atoms,
# pyrantel 14; meglumine has more atoms than diatrizoic acid once hydrogens count),
# so "keep the largest fragment" alone returns the counter-ion. Each is passed
# over while the input has another organic component. Not listed, on purpose:
# partners that are drugs in their own right (salicylic and benzoic acid,
# piperazine, procaine); leaving one of those out is reported instead.
_PARTNERS = frozenset(
    _flat(Chem.MolFromSmiles(smiles))
    for smiles in (
        # acids
        "OC(=O)c1cc2ccccc2c(Cc2c(O)c(cc3ccccc23)C(O)=O)c1O",  # pamoic (embonic)
        "OC(C(O)C(O)=O)C(O)=O",  # tartaric
        "OC(=O)CC(O)(CC(O)=O)C(O)=O",  # citric
        "OC(=O)C=CC(O)=O",  # maleic, fumaric
        "OC(=O)CCC(O)=O",  # succinic
        "OC(CC(O)=O)C(O)=O",  # malic
        "CC(O)C(O)=O",  # lactic
        "OCC(O)C(O)C(O)C(O)C(O)=O",  # gluconic
        "OC(C(O)C(O)C(O)=O)C(O)C(O)=O",  # glucaric (saccharic), mucic
        "OC(=O)CCCCC(O)=O",  # adipic
        "OC(=O)CNC(=O)c1ccccc1",  # hippuric
        "OC(=O)C(O)=O",  # oxalic
        "OC(=O)CC(O)=O",  # malonic
        "CC(O)=O",  # acetic
        "OC=O",  # formic
        "OC(=O)C(F)(F)F",  # trifluoroacetic
        "CS(O)(=O)=O",  # methanesulfonic (mesylate)
        "CCS(O)(=O)=O",  # ethanesulfonic (esylate)
        "OS(=O)(=O)CCS(O)(=O)=O",  # ethanedisulfonic (edisylate)
        "OCCS(O)(=O)=O",  # isethionic
        "OS(=O)(=O)c1ccccc1",  # benzenesulfonic (besylate)
        "Cc1ccc(cc1)S(O)(=O)=O",  # p-toluenesulfonic (tosylate)
        "OS(=O)(=O)c1ccc2ccccc2c1",  # naphthalene-2-sulfonic (napsylate)
        "OS(=O)(=O)c1cccc2c1cccc2S(O)(=O)=O",  # naphthalene-1,5-disulfonic
        "CC1(C)C2CCC1(CS(O)(=O)=O)C(=O)C2",  # camphorsulfonic (camsylate)
        "OC(=O)c1ccc2ccccc2c1O",  # 1-hydroxy-2-naphthoic (xinafoate)
        "CCCCCCCCCCCCCCCCCC(O)=O",  # stearic
        "CCCCCCCCCCCCOS(O)(=O)=O",  # lauryl sulfate (estolate)
        # bases
        "CNCC(O)C(O)C(O)C(O)CO",  # meglumine
        "NC(CO)(CO)CO",  # tromethamine
        "C[N+](C)(C)CCO",  # choline
        "OCCN(CCO)CCO",  # trolamine
        "OCCNCCO",  # diolamine
        "NCCO",  # olamine
        "C(CNCc1ccccc1)NCc1ccccc1",  # benzathine
        "NCCCCC(N)C(O)=O",  # lysine
        "NC(CCCNC(N)=N)C(O)=O",  # arginine
        "CCNCC",  # diethylamine
        "CC(C)(C)N",  # tert-butylamine (erbumine)
        "C1CCC(CC1)NC1CCCCC1",  # dicyclohexylamine
        "OCCN1CCCC1",  # epolamine
        "NCCN",  # ethylenediamine
        # solvents
        "CO",
        "CCO",
        "CC(C)O",
        "CC(O)CO",  # propylene glycol
        "CC(C)=O",
        "CS(C)=O",
        "CC#N",
        "ClCCl",
        "ClC(Cl)Cl",
        "CCOC(C)=O",
        "CCOCC",
        "C1CCOC1",
        "C1COCCO1",
        "CN(C)C=O",
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
    2. Keep the organic fragment with the most heavy atoms, passing over common
       salt formers and solvents (pamoate, meglumine, tartrate, ethanol) while
       there is another. Everything else is reported as removed, with a note
       when it includes an organic component that is not one of those.
    3. Neutralize charges (optional, default on).

    Two kinds of input have no parent to extract and are returned as given, with
    a note: compounds of a metal other than groups 1 and 2 (cisplatin, ferrocene,
    ferrous fumarate: without the metal a ligand or a counter-ion is left), and
    inputs without an organic part (sodium chloride).

    Example:
        >>> MoleculeStandardizer().standardize_smiles("CC(=O)[O-].[Na+]")
        'CC(=O)O'
    """

    def __init__(self, neutralize: bool = True) -> None:
        self.neutralize = neutralize
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
        metals = _compound_metals(mol)
        if metals:
            return Standardization(
                Chem.Mol(mol),
                note=f"Contains {_listed(metals)}: without the metal only a ligand "
                "or a counter-ion would be left, so the input is returned as given.",
            )
        clean = _neutral_sulfoxides(rdMolStandardize.Cleanup(mol))
        fragments = list(Chem.GetMolFrags(clean, asMols=True))
        organic = [frag for frag in fragments if _has_carbon(frag)]
        if not organic:
            return Standardization(
                Chem.Mol(mol),
                note="No organic component, so there is no parent to extract: "
                "returned as given.",
            )
        candidates = [frag for frag in organic if _flat(frag) not in _PARTNERS]
        chosen = max(candidates or organic, key=_size)
        parent = self._uncharger.uncharge(chosen) if self.neutralize else chosen
        Chem.SanitizeMol(parent)
        # A 2:1 salt holds the parent twice: the second copy was not removed.
        same = _flat(chosen)
        others = [
            frag for frag in fragments if frag is not chosen and _flat(frag) != same
        ]
        return Standardization(
            parent,
            removed=[str(Chem.MolToSmiles(frag)) for frag in others],
            note=_note(others, only_partners=not candidates and _has_carbon_in(others)),
        )

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


def _has_carbon_in(fragments: list[Chem.Mol]) -> bool:
    return any(_has_carbon(frag) for frag in fragments)


def _listed(items: list[str]) -> str:
    return ", ".join(items)


def _size(mol: Chem.Mol) -> tuple[int, float, str]:
    """Heavy atoms first: counting hydrogens favours sugar-like counter-ions."""
    weight = float(rdMolDescriptors._CalcMolWt(mol))
    return mol.GetNumHeavyAtoms(), weight, str(Chem.MolToSmiles(mol))


def _compound_metals(mol: Chem.Mol) -> list[str]:
    """Metals that make ``mol`` a metal compound rather than a salt: any outside
    groups 1 and 2, drawn bonded or as an ion (PubChem writes carboplatin as
    ``[NH2-].[NH2-].[Pt+2]`` and a ligand), and any metal bonded to carbon."""
    found = {
        atom.GetSymbol()
        for atom in mol.GetAtoms()
        if _is_metal(atom)
        and (
            atom.GetAtomicNum() not in _SALT_METALS
            or any(other.GetAtomicNum() == 6 for other in atom.GetNeighbors())
        )
    }
    return sorted(found)


def _note(removed: list[Chem.Mol], only_partners: bool) -> str | None:
    """What a reader has to know about the fragments that were left out."""
    notes = []
    metals = sorted(
        {
            atom.GetSymbol()
            for frag in removed
            for atom in frag.GetAtoms()
            if _is_metal(atom) and atom.GetAtomicNum() not in _QUIET_METALS
        }
    )
    if metals:
        notes.append(
            f"Left out {_listed(metals)}: if the metal is the active part, use the "
            "input as it is."
        )
    unknown = sorted(
        {
            str(Chem.MolToSmiles(frag))
            for frag in removed
            if _has_carbon(frag) and _flat(frag) not in _PARTNERS
        }
    )
    if unknown:
        notes.append(
            f"Left out {_listed(unknown)}, which is not a common counter-ion or "
            "solvent: the input may be a mixture, and the parent returned is only "
            "its largest component."
        )
    if only_partners:
        notes.append(
            "Every organic component is a common counter-ion or solvent; the "
            "largest was kept."
        )
    return " ".join(notes) or None


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
