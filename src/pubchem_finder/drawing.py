"""Render molecules to PNG with RDKit, optionally highlighting what matters.

Highlight modes:

- ``smarts``: atoms and bonds matching a SMARTS pattern (e.g. a pharmacophore).
- ``diff``: for each molecule, the atoms that are *not* part of its maximum
  common substructure with a reference: "what changed in this analogue".
  Depictions are aligned to the shared core so the changes line up visually.
  Rings are shared whole or not at all, and an aromatic ring is not a saturated
  one. A molecule that shares nothing with the reference is highlighted whole. A
  panel whose common substructure could not be found in time is drawn without
  highlights and says so in its caption.
- ``similarity_map``: atom contributions to ECFP4 similarity with a reference
  (green raises similarity, pink lowers it).
"""

import io
from collections.abc import Sequence
from typing import Literal

from PIL import Image as PILImage
from PIL import ImageDraw, ImageFont
from rdkit import Chem, DataStructs
from rdkit.Chem import (
    Draw,
    rdDepictor,
    rdFingerprintGenerator,
    rdFMCS,
    rdRascalMCES,
)
from rdkit.Chem.Draw import SimilarityMaps, rdMolDraw2D

from .core.exceptions import InvalidInput, MoleculeValidationError
from .core.limits import check_atoms, check_length

Highlight = Literal["none", "smarts", "diff", "similarity_map"]
MAX_MOLECULES = 24
MCS_TIMEOUT_S = 2
NOT_COMPUTED = "difference not computed"
MACROCYCLE_ATOMS = 12
# The fallback search, when FindMCS runs out of time (see _shared_atoms).
_RASCAL_OPTIONS: dict[str, object] = {
    "similarityThreshold": 0.0,  # the caller chose the pair; compare any two
    "ringMatchesRingOnly": True,
    "completeAromaticRings": True,
    "completeSmallestRings": True,  # as FindMCS's completeRingsOnly
    "singleLargestFrag": True,
    "timeout": MCS_TIMEOUT_S,
}
# RDKit's option objects take any attribute name without complaint, and RDKit
# 2024.09 has no completeSmallestRings: RASCAL then shares part of a ring, which
# FindMCS never does, so with such an RDKit the fallback is not used at all.
_RASCAL_USABLE = all(
    hasattr(rdRascalMCES.RascalOptions, name) for name in _RASCAL_OPTIONS
)
# Laying a long chain out along the shared core makes its similarity map cost grow
# with the cube of its size while holding the GIL (18 s at 456 atoms, 0.8 s with
# RDKit's own layout); larger molecules keep their own layout.
ALIGN_MAX_ATOMS = 100
# Even unaligned, a map costs about the square of a molecule's size (5 s at 1,000
# atoms), so one similarity-map request draws at most this many atoms in total.
MAX_SIMILARITY_MAP_ATOMS = 2_500


def _parse(smiles: Sequence[str]) -> list[Chem.Mol]:
    mols: list[Chem.Mol] = []
    for i, s in enumerate(smiles):
        check_length(s.strip(), f"SMILES at position {i + 1}")
        mol = Chem.MolFromSmiles(s.strip()) if s.strip() else None
        if mol is None:
            raise MoleculeValidationError(f"Invalid SMILES at position {i + 1}: {s!r}")
        mols.append(check_atoms(mol, f"The molecule at position {i + 1}"))
    return mols


def _bonds_between(mol: Chem.Mol, atoms: set[int]) -> list[int]:
    return [
        b.GetIdx()
        for b in mol.GetBonds()
        if b.GetBeginAtomIdx() in atoms and b.GetEndAtomIdx() in atoms
    ]


def _bonds_touching(mol: Chem.Mol, atoms: set[int]) -> list[int]:
    return [
        b.GetIdx()
        for b in mol.GetBonds()
        if b.GetBeginAtomIdx() in atoms or b.GetEndAtomIdx() in atoms
    ]


def smarts_highlights(
    mols: Sequence[Chem.Mol], smarts: str
) -> tuple[list[list[int]], list[list[int]]]:
    """Atoms and bonds of every match of ``smarts`` in each molecule."""
    check_length(smarts, "SMARTS")
    pattern = Chem.MolFromSmarts(smarts)
    if pattern is None:
        raise InvalidInput(f"Invalid SMARTS: {smarts!r}")
    check_atoms(pattern, "SMARTS")
    atom_lists: list[list[int]] = []
    bond_lists: list[list[int]] = []
    for mol in mols:
        atoms = {a for match in mol.GetSubstructMatches(pattern) for a in match}
        atom_lists.append(sorted(atoms))
        bond_lists.append(_bonds_between(mol, atoms))
    return atom_lists, bond_lists


def diff_highlights(
    mols: Sequence[Chem.Mol], reference: Chem.Mol
) -> tuple[list[list[int]], list[list[int]], list[int]]:
    """Atoms outside each molecule's MCS with ``reference``, and their bonds.

    Also aligns each depiction to the shared core of the reference. A molecule
    that shares nothing with the reference has every atom highlighted. The third
    list holds the positions of molecules whose common substructure could not be
    found in time: they get no highlights, which must not be read as "the same".
    """
    rdDepictor.Compute2DCoords(reference)
    atom_lists: list[list[int]] = []
    bond_lists: list[list[int]] = []
    not_computed: list[int] = []
    for position, mol in enumerate(mols):
        shared = _shared_atoms(reference, mol)
        if shared is None:
            not_computed.append(position)
            atom_lists.append([])
            bond_lists.append([])
            continue
        if shared:
            try:
                rdDepictor.GenerateDepictionMatching2DStructure(
                    mol, reference, atomMap=shared
                )
            except ValueError:
                rdDepictor.Compute2DCoords(mol)
        changed = {a.GetIdx() for a in mol.GetAtoms()} - {m for _, m in shared}
        atom_lists.append(sorted(changed))
        bond_lists.append(_bonds_touching(mol, changed))
    return atom_lists, bond_lists, not_computed


def _shared_atoms(reference: Chem.Mol, mol: Chem.Mol) -> list[tuple[int, int]] | None:
    """(reference atom, molecule atom) pairs of the maximum common substructure;
    empty when the two share none, None when it was not found in time.

    Bond orders have to match exactly: by FindMCS's default an aromatic bond
    matches a single or a double one, and a phenyl ring that became cyclohexyl
    counted as unchanged.

    FindMCS stops at its time limit with whatever it has: for fused-ring analogues
    such as dolutegravir and cabotegravir that is 2 atoms of a 25-atom core, and
    about one drug-like analogue pair in eleven runs into the limit. A search
    that was cut short is therefore never used; RASCAL is asked instead.
    """
    result = rdFMCS.FindMCS(
        [reference, mol],
        timeout=MCS_TIMEOUT_S,
        ringMatchesRingOnly=True,
        completeRingsOnly=True,
        bondCompare=rdFMCS.BondCompare.CompareOrderExact,
    )
    if result.canceled:
        return _rascal_atoms(reference, mol, at_least=result.numAtoms)
    core = Chem.MolFromSmarts(result.smartsString) if result.numAtoms else None
    if core is None:
        return []
    ref_match = reference.GetSubstructMatch(core)
    mol_match = mol.GetSubstructMatch(core)
    if not ref_match or not mol_match:
        return []
    return list(zip(ref_match, mol_match, strict=True))


def _rascal_atoms(
    reference: Chem.Mol, mol: Chem.Mol, at_least: int
) -> list[tuple[int, int]] | None:
    """The common substructure by RASCAL, or None when its answer cannot be relied
    on. ``at_least`` is the size FindMCS had reached when it was stopped.

    For close analogues without large rings RASCAL finds FindMCS's core in
    milliseconds (of 21 analogue pairs FindMCS needed more than its limit for, 17
    gave the same atoms and 2 differed by one). It also gives up without saying
    so: on a macrocycle it returned 12 atoms of a 65-atom core, not flagged as
    timed out, after several seconds. So no macrocycles, and no answer smaller
    than FindMCS's partial one or than half of the smaller molecule.
    """
    if not _RASCAL_USABLE or _has_macrocycle(reference) or _has_macrocycle(mol):
        return None
    # RDKit's stubs type every option as RascalOptions, hence setattr.
    options = rdRascalMCES.RascalOptions()
    for name, value in _RASCAL_OPTIONS.items():
        setattr(options, name, value)
    found = rdRascalMCES.FindMCES(reference, mol, options)
    if not found or found[0].timedOut:
        return None
    pairs = [(int(r), int(m)) for r, m in found[0].atomMatches()]
    smaller = min(reference.GetNumAtoms(), mol.GetNumAtoms())
    if len(pairs) < at_least or 2 * len(pairs) < smaller:
        return None
    return pairs


def _has_macrocycle(mol: Chem.Mol) -> bool:
    return any(len(ring) >= MACROCYCLE_ATOMS for ring in mol.GetRingInfo().AtomRings())


_MORGAN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def atom_contributions(reference: Chem.Mol, probe: Chem.Mol) -> list[float]:
    """How much each probe atom adds to its ECFP4 Tanimoto similarity with
    ``reference``: similarity minus similarity without the bits centred on it."""
    ref_fp = _MORGAN.GetFingerprint(reference)
    atoms = list(range(probe.GetNumAtoms()))
    full = DataStructs.TanimotoSimilarity(ref_fp, _MORGAN.GetFingerprint(probe))
    return [
        full
        - DataStructs.TanimotoSimilarity(
            ref_fp,
            _MORGAN.GetFingerprint(probe, fromAtoms=[a for a in atoms if a != i]),
        )
        for i in atoms
    ]


LEGEND_HEIGHT = 28
# RDKit shrinks atom labels along with the molecule, down to 6 px by default,
# which makes drug-sized molecules unreadable in a 300 px panel.
MIN_FONT_SIZE = 13


def _draw_options() -> Draw.MolDrawOptions:
    options = Draw.MolDrawOptions()
    options.minFontSize = MIN_FONT_SIZE  # type: ignore[assignment]  # stub says MolDrawOptions
    return options


def _similarity_map_png(
    reference: Chem.Mol, probe: Chem.Mol, size: int, legend: str
) -> bytes:
    map_height = size - LEGEND_HEIGHT if legend else size
    drawer = rdMolDraw2D.MolDraw2DCairo(size, map_height)
    drawer.drawOptions().minFontSize = MIN_FONT_SIZE  # type: ignore[assignment]
    weights = atom_contributions(reference, probe)
    if any(weights):
        weights, _ = SimilarityMaps.GetStandardizedWeights(weights)  # type: ignore[no-untyped-call]
    SimilarityMaps.GetSimilarityMapFromWeights(probe, weights, draw2d=drawer)  # type: ignore[no-untyped-call]
    drawer.FinishDrawing()
    panel = PILImage.new("RGB", (size, size), "white")
    panel.paste(
        PILImage.open(io.BytesIO(drawer.GetDrawingText())).convert("RGB"), (0, 0)
    )
    if legend:
        # Pillow draws captions at a fixed size; RDKit would scale them with the molecule.
        font = ImageFont.load_default(size=16)
        ImageDraw.Draw(panel).text(
            (size / 2, size - LEGEND_HEIGHT / 2),
            legend,
            fill="black",
            font=font,
            anchor="mm",
        )
    out = io.BytesIO()
    panel.save(out, format="PNG")
    return out.getvalue()


def _tile(pngs: Sequence[bytes], per_row: int, size: int) -> bytes:
    rows = (len(pngs) + per_row - 1) // per_row
    cols = min(per_row, len(pngs))
    canvas = PILImage.new("RGB", (cols * size, rows * size), "white")
    for i, data in enumerate(pngs):
        tile = PILImage.open(io.BytesIO(data)).convert("RGB")
        canvas.paste(tile, ((i % per_row) * size, (i // per_row) * size))
    out = io.BytesIO()
    canvas.save(out, format="PNG", optimize=True)
    return out.getvalue()


def draw_grid(
    smiles: Sequence[str],
    legends: Sequence[str] | None = None,
    highlight: Highlight = "none",
    smarts: str | None = None,
    reference: str | None = None,
    size: int = 300,
    per_row: int = 4,
) -> bytes:
    """
    Draw up to 24 molecules as one PNG grid.

    Args:
        smiles: Molecules to draw.
        legends: One caption per molecule (defaults to none).
        highlight: ``none``, ``smarts`` (needs ``smarts``), ``diff`` or
            ``similarity_map`` (both need ``reference``).
        smarts: Pattern to highlight in ``smarts`` mode.
        reference: SMILES to compare against in ``diff`` / ``similarity_map``.
        size: Pixel width and height of each molecule panel (150-600).
        per_row: Panels per row (1-8).

    Returns:
        PNG bytes.
    """
    if not 1 <= len(smiles) <= MAX_MOLECULES:
        raise InvalidInput(f"Draw between 1 and {MAX_MOLECULES} molecules")
    if not 150 <= size <= 600:
        raise InvalidInput("size must be between 150 and 600 pixels")
    if not 1 <= per_row <= 8:
        raise InvalidInput("per_row must be between 1 and 8")
    if legends is not None and len(legends) != len(smiles):
        raise InvalidInput("Give one legend per molecule, or none")
    mols = _parse(smiles)
    labels = list(legends) if legends is not None else [""] * len(mols)

    if highlight in ("diff", "similarity_map"):
        if not reference:
            raise InvalidInput(f"highlight={highlight!r} needs a reference SMILES")
        (ref,) = _parse([reference])
    if highlight == "smarts" and not smarts:
        raise InvalidInput("highlight='smarts' needs a SMARTS pattern")

    if highlight == "similarity_map":
        total = sum(m.GetNumAtoms() for m in mols)
        if total > MAX_SIMILARITY_MAP_ATOMS:
            raise InvalidInput(
                f"Similarity maps of {total:,} atoms in one request; the limit is "
                f"{MAX_SIMILARITY_MAP_ATOMS:,}. Draw fewer or smaller molecules."
            )
        # Only to align depictions to the shared core.
        diff_highlights([m for m in mols if m.GetNumAtoms() <= ALIGN_MAX_ATOMS], ref)
        panels = [
            _similarity_map_png(ref, mol, size, label)
            for mol, label in zip(mols, labels, strict=True)
        ]
        return _tile(panels, per_row, size)

    atom_lists: list[list[int]] | None = None
    bond_lists: list[list[int]] | None = None
    if highlight == "smarts":
        atom_lists, bond_lists = smarts_highlights(mols, smarts or "")
    elif highlight == "diff":
        atom_lists, bond_lists, not_computed = diff_highlights(mols, ref)
        for position in not_computed:
            label = labels[position]
            labels[position] = f"{label} ({NOT_COMPUTED})" if label else NOT_COMPUTED
    elif highlight != "none":
        raise InvalidInput(f"Unknown highlight mode {highlight!r}")

    png = Draw.MolsToGridImage(  # type: ignore[no-untyped-call]
        mols,
        molsPerRow=min(per_row, len(mols)),
        subImgSize=(size, size),
        legends=labels,
        highlightAtomLists=atom_lists,
        highlightBondLists=bond_lists,
        returnPNG=True,
        drawOptions=_draw_options(),
    )
    return bytes(png)
