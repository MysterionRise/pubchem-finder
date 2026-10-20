"""Molecule standardization."""

import logging

from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

logger = logging.getLogger(__name__)


class MoleculeStandardizer:
    """
    Standardize molecules to a canonical form before fingerprinting.

    Pipeline (the order matters: the parent fragment is chosen before
    neutralization so counter-ions never influence the charge balance):

    1. ``rdMolStandardize.Cleanup``: remove explicit hydrogens, normalize
       functional groups, reionize, sanitize.
    2. Keep the largest organic fragment (drops salts and counter-ions).
    3. Neutralize charges (optional, default on).

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
        if mol is None:
            return None
        try:
            mol = rdMolStandardize.Cleanup(mol)
            mol = self._fragment_chooser.choose(mol)
            if self.neutralize:
                mol = self._uncharger.uncharge(mol)
            Chem.SanitizeMol(mol)
        except Exception:
            logger.warning("Standardization failed", exc_info=True)
            return None
        return mol

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
