"""Core molecule data structures."""

from dataclasses import dataclass, field
from datetime import datetime
from typing import Any, Dict, Optional

from .exceptions import MoleculeValidationError


@dataclass
class Molecule:
    """
    Core molecule data structure.

    This is the central data model used throughout the application.
    Intentionally simple and focused on essentials.

    Attributes:
        id: Unique identifier (UUID or PubChem CID)
        smiles: Original SMILES string
        canonical_smiles: Canonicalized SMILES string
        inchi: InChI string
        inchi_key: InChI Key (27 characters)
        molecular_weight: Molecular weight in Daltons
        logp: Calculated LogP (octanol-water partition coefficient)
        num_h_donors: Number of hydrogen bond donors
        num_h_acceptors: Number of hydrogen bond acceptors
        tpsa: Topological polar surface area
        num_rotatable_bonds: Number of rotatable bonds
        num_aromatic_rings: Number of aromatic rings
        fingerprints: Dictionary of fingerprints {type: vector}
        source: Data source (e.g., "pubchem", "chembl")
        metadata: Additional metadata
        created_at: Creation timestamp
        updated_at: Last update timestamp
    """

    # Identifiers
    id: Optional[str] = None

    # Chemical representations
    smiles: str = ""
    canonical_smiles: str = ""
    inchi: Optional[str] = None
    inchi_key: Optional[str] = None

    # Molecular properties
    molecular_weight: Optional[float] = None
    logp: Optional[float] = None
    num_h_donors: Optional[int] = None
    num_h_acceptors: Optional[int] = None
    tpsa: Optional[float] = None
    num_rotatable_bonds: Optional[int] = None
    num_aromatic_rings: Optional[int] = None

    # Fingerprints (stored as dict: {type: vector})
    fingerprints: Dict[str, Any] = field(default_factory=dict)

    # Metadata
    source: str = "pubchem"
    metadata: Dict[str, Any] = field(default_factory=dict)
    created_at: Optional[datetime] = None
    updated_at: Optional[datetime] = None

    def __post_init__(self) -> None:
        """Validate essential fields."""
        if not self.smiles and not self.canonical_smiles:
            raise MoleculeValidationError(
                "Either smiles or canonical_smiles must be provided"
            )

        if self.created_at is None:
            self.created_at = datetime.utcnow()

    def add_fingerprint(self, fp_type: str, fingerprint: Any) -> None:
        """
        Add a fingerprint to the molecule.

        Args:
            fp_type: Fingerprint type identifier
            fingerprint: Fingerprint data (numpy array or other format)
        """
        self.fingerprints[fp_type] = fingerprint

    def get_fingerprint(self, fp_type: str) -> Optional[Any]:
        """
        Get a fingerprint by type.

        Args:
            fp_type: Fingerprint type identifier

        Returns:
            Fingerprint data or None if not found
        """
        return self.fingerprints.get(fp_type)

    def has_fingerprint(self, fp_type: str) -> bool:
        """Check if molecule has a specific fingerprint type."""
        return fp_type in self.fingerprints

    def __repr__(self) -> str:
        """String representation."""
        return (
            f"Molecule(id={self.id}, "
            f"smiles={self.canonical_smiles or self.smiles}, "
            f"mw={self.molecular_weight:.2f if self.molecular_weight else 'N/A'})"
        )
