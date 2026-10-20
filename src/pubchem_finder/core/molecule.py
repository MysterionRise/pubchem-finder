"""Core molecule data structure."""

from dataclasses import dataclass, field
from datetime import UTC, datetime
from typing import Any

from .exceptions import MoleculeValidationError


@dataclass(eq=False)
class Molecule:
    """
    Central molecule record used throughout the application.

    Equality is identity-based (``eq=False``) because ``fingerprints`` holds
    numpy arrays, whose ``==`` does not reduce to a single boolean.

    Attributes:
        id: Unique identifier (for PubChem records, the CID as a string)
        smiles: Original SMILES string
        canonical_smiles: Canonical SMILES string
        inchi: InChI string
        inchi_key: InChIKey (27 characters)
        molecular_weight: Molecular weight in Daltons
        logp: Crippen LogP
        num_h_donors: Number of hydrogen bond donors
        num_h_acceptors: Number of hydrogen bond acceptors
        tpsa: Topological polar surface area
        num_rotatable_bonds: Number of rotatable bonds
        num_aromatic_rings: Number of aromatic rings
        num_heavy_atoms: Number of heavy (non-hydrogen) atoms
        fingerprints: Fingerprints keyed by type name
        source: Data source (e.g. ``"pubchem"``)
        metadata: Free-form extra data
        created_at: Creation timestamp (UTC)
        updated_at: Last update timestamp (UTC)
    """

    # Identifiers
    id: str | None = None

    # Chemical representations
    smiles: str = ""
    canonical_smiles: str = ""
    inchi: str | None = None
    inchi_key: str | None = None

    # Molecular properties
    molecular_weight: float | None = None
    logp: float | None = None
    num_h_donors: int | None = None
    num_h_acceptors: int | None = None
    tpsa: float | None = None
    num_rotatable_bonds: int | None = None
    num_aromatic_rings: int | None = None
    num_heavy_atoms: int | None = None

    # Fingerprints keyed by type name
    fingerprints: dict[str, Any] = field(default_factory=dict)

    # Metadata
    source: str = "pubchem"
    metadata: dict[str, Any] = field(default_factory=dict)
    created_at: datetime | None = None
    updated_at: datetime | None = None

    def __post_init__(self) -> None:
        if not self.smiles and not self.canonical_smiles:
            raise MoleculeValidationError(
                "Either smiles or canonical_smiles must be provided"
            )
        if self.created_at is None:
            self.created_at = datetime.now(UTC)

    def add_fingerprint(self, fp_type: str, fingerprint: Any) -> None:
        """Store a fingerprint under ``fp_type``."""
        self.fingerprints[fp_type] = fingerprint

    def get_fingerprint(self, fp_type: str) -> Any | None:
        """Return the fingerprint stored under ``fp_type``, or ``None``."""
        return self.fingerprints.get(fp_type)

    def has_fingerprint(self, fp_type: str) -> bool:
        """Whether a fingerprint of type ``fp_type`` is stored."""
        return fp_type in self.fingerprints

    def __repr__(self) -> str:
        mw = "N/A" if self.molecular_weight is None else f"{self.molecular_weight:.2f}"
        return (
            f"Molecule(id={self.id}, "
            f"smiles={self.canonical_smiles or self.smiles}, "
            f"mw={mw})"
        )
