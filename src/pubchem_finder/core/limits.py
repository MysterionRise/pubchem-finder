"""Size limits for structures that users and MCP clients pass in.

RDKit's SMILES writer is recursive and its cost grows with the square of a chain's
length: an 80,000-character chain overflows a worker thread's stack and takes the
MCP server down with it, and long chains make Murcko scaffolds and similarity maps
run for minutes while holding the GIL, which stalls an HTTP server for every client.
The limits sit well above real compounds: the largest drug in the demo set has 492
heavy atoms and a 1,185-character SMILES.
"""

import re

from rdkit import Chem

from .exceptions import MoleculeValidationError

MAX_STRUCTURE_CHARS = 10_000
MAX_ATOMS = 1_000

_CX_EXTENSION = re.compile(r"\s\|")
_ATOMIC_NUMBER_ATOM = re.compile(r"\[#\d")


def check_length(text: str, what: str = "SMILES") -> None:
    """Refuse structure text longer than :data:`MAX_STRUCTURE_CHARS`, before parsing."""
    if len(text) > MAX_STRUCTURE_CHARS:
        raise MoleculeValidationError(
            f"{what} is {len(text):,} characters long; the limit is "
            f"{MAX_STRUCTURE_CHARS:,}"
        )


def check_one_structure(text: str, what: str = "SMILES") -> None:
    """Refuse text that is not a single structure, before parsing.

    RDKit reads ``SMILES name`` as a molecule name and stops at the first non-ASCII
    character, so anything after a space, tab, line break or such a character would
    be dropped without a word. Pass the stripped text.

    A CXSMILES extension (``C[C@H](N)C(=O)O |&1:1|``) is such text too: RDKit
    reads it, but nothing here can write it back, so a racemate came back as one
    enantiomer. It gets a message of its own.
    """
    if _CX_EXTENSION.search(text):
        raise MoleculeValidationError(
            f"{what} has a CXSMILES extension (the |...| part after the space), which "
            "is not supported: enhanced stereo such as |&1:1|, polymer brackets and "
            "attachment points would be dropped. Send the structure without it."
        )
    if any(not "!" <= c <= "~" for c in text):
        raise MoleculeValidationError(
            f"{what} must be one structure with no spaces, line breaks or non-ASCII "
            "characters; send one structure per call and join the parts of a mixture "
            "with '.'"
        )


def check_smiles_atoms(text: str, what: str = "SMILES") -> None:
    """Refuse SMARTS atoms in text that is about to be read as SMILES.

    RDKit reads ``[#6]`` in a SMILES as a carbon without hydrogens, so a SMARTS
    passed where a SMILES belongs came back as a radical: phenol written
    ``[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#8]`` weighed 88.07 instead of 94.11,
    with no warning.
    """
    if _ATOMIC_NUMBER_ATOM.search(text):
        raise MoleculeValidationError(
            f"{what} has atoms written as atomic numbers (such as [#6]), which is "
            "SMARTS: read as SMILES they would be atoms without hydrogens. Write "
            "element symbols (c1ccccc1O), or search with substructure_search and "
            'query_type="smarts" (`substructure --smarts` on the command line).'
        )


def check_atoms(mol: Chem.Mol, what: str = "Molecule") -> Chem.Mol:
    """Refuse a parsed molecule or query with more than :data:`MAX_ATOMS` atoms."""
    n = mol.GetNumAtoms()
    if n > MAX_ATOMS:
        raise MoleculeValidationError(
            f"{what} has {n:,} atoms; the limit is {MAX_ATOMS:,}"
        )
    return mol
