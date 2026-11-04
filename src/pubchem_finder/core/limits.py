"""Size limits for structures that users and MCP clients pass in.

RDKit's SMILES writer is recursive and its cost grows with the square of a chain's
length: an 80,000-character chain overflows a worker thread's stack and takes the
MCP server down with it, and long chains make Murcko scaffolds and similarity maps
run for minutes while holding the GIL, which stalls an HTTP server for every client.
The limits sit well above real compounds: the largest drug in the demo set has 492
heavy atoms and a 1,185-character SMILES.
"""

from rdkit import Chem

from .exceptions import MoleculeValidationError

MAX_STRUCTURE_CHARS = 10_000
MAX_ATOMS = 1_000


def check_length(text: str, what: str = "SMILES") -> None:
    """Refuse structure text longer than :data:`MAX_STRUCTURE_CHARS`, before parsing."""
    if len(text) > MAX_STRUCTURE_CHARS:
        raise MoleculeValidationError(
            f"{what} is {len(text):,} characters long; the limit is "
            f"{MAX_STRUCTURE_CHARS:,}"
        )


def check_atoms(mol: Chem.Mol, what: str = "Molecule") -> Chem.Mol:
    """Refuse a parsed molecule or query with more than :data:`MAX_ATOMS` atoms."""
    n = mol.GetNumAtoms()
    if n > MAX_ATOMS:
        raise MoleculeValidationError(
            f"{what} has {n:,} atoms; the limit is {MAX_ATOMS:,}"
        )
    return mol
