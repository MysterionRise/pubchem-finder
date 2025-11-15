"""Molecule standardization and normalization."""

import logging
from typing import Optional

from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

logger = logging.getLogger(__name__)


class MoleculeStandardizer:
    """
    Standardize molecules to canonical form.

    This is CRITICAL for consistent fingerprints and searches.
    All molecules should be standardized before fingerprinting.

    The standardization process:
    1. Remove hydrogen atoms
    2. Normalize functional groups
    3. Neutralize charges (optional)
    4. Select largest fragment (remove salts/counterions)
    5. Canonicalize

    Example:
        >>> standardizer = MoleculeStandardizer()
        >>> mol = Chem.MolFromSmiles("CCO.Cl")  # Ethanol with chloride salt
        >>> std_mol = standardizer.standardize(mol)
        >>> Chem.MolToSmiles(std_mol)
        'CCO'
    """

    def __init__(self, neutralize: bool = True):
        """
        Initialize standardizer.

        Args:
            neutralize: Whether to neutralize charges (default: True)
        """
        self.neutralize = neutralize
        self.normalizer = rdMolStandardize.Normalizer()
        self.uncharger = rdMolStandardize.Uncharger()

    def standardize(self, mol: Chem.Mol) -> Optional[Chem.Mol]:
        """
        Standardize a molecule.

        Args:
            mol: RDKit Mol object

        Returns:
            Standardized molecule or None if standardization fails
        """
        if mol is None:
            return None

        try:
            # Remove explicit hydrogens
            mol = Chem.RemoveHs(mol)

            # Normalize (standardize functional groups)
            mol = self.normalizer.normalize(mol)

            # Remove charges (makes similarity more robust)
            if self.neutralize:
                mol = self.uncharger.uncharge(mol)

            # Take largest fragment (removes salts, counterions)
            mol = self._get_largest_fragment(mol)

            # Canonicalize
            Chem.SanitizeMol(mol)

            return mol

        except Exception as e:
            logger.error(f"Standardization failed: {e}")
            return None

    def _get_largest_fragment(self, mol: Chem.Mol) -> Chem.Mol:
        """
        Select largest fragment from molecule.

        If molecule contains multiple fragments (e.g., salt),
        returns the fragment with the most heavy atoms.

        Args:
            mol: RDKit Mol object

        Returns:
            Largest fragment
        """
        frags = Chem.GetMolFrags(mol, asMols=True)
        if len(frags) == 1:
            return mol
        # Return fragment with most heavy atoms
        return max(frags, key=lambda m: m.GetNumHeavyAtoms())

    def standardize_smiles(self, smiles: str) -> Optional[str]:
        """
        Standardize SMILES string.

        Args:
            smiles: SMILES string

        Returns:
            Standardized SMILES or None if standardization fails
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.error(f"Invalid SMILES: {smiles}")
            return None

        mol = self.standardize(mol)
        if mol is None:
            return None

        return Chem.MolToSmiles(mol)
