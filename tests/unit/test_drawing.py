"""Tests for molecule drawing and highlighting."""

import io

import pytest
from PIL import Image
from rdkit import Chem

from pubchem_finder import drawing
from pubchem_finder.core.exceptions import InvalidInput, MoleculeValidationError

ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"
SALICYLIC = "O=C(O)c1ccccc1O"
METHYL_SALICYLATE = "COC(=O)c1ccccc1O"
PNG_MAGIC = b"\x89PNG\r\n\x1a\n"


def image(png: bytes) -> Image.Image:
    assert png.startswith(PNG_MAGIC)
    return Image.open(io.BytesIO(png))


def test_plain_grid_size():
    img = image(drawing.draw_grid([ASPIRIN, SALICYLIC, "CCO"], size=200, per_row=2))
    assert img.size == (400, 400)  # 2 columns x 2 rows of 200 px


def test_single_row_uses_actual_count():
    img = image(drawing.draw_grid([ASPIRIN, SALICYLIC], size=200, per_row=4))
    assert img.size == (400, 200)


@pytest.mark.parametrize(
    ("mode", "kwargs"),
    [
        ("smarts", {"smarts": "c1ccccc1C(=O)[OX2H1]"}),
        ("diff", {"reference": ASPIRIN}),
        ("similarity_map", {"reference": ASPIRIN}),
    ],
)
def test_highlight_modes_render(mode, kwargs):
    png = drawing.draw_grid(
        [ASPIRIN, METHYL_SALICYLATE],
        ["aspirin", "methyl salicylate"],
        highlight=mode,
        **kwargs,
    )
    assert image(png).size[0] > 0


def test_smarts_highlights_match_atoms():
    mols = [Chem.MolFromSmiles(s) for s in (SALICYLIC, "CCO")]
    atoms, bonds = drawing.smarts_highlights(mols, "[OX2H]")
    assert len(atoms[0]) == 2  # the phenol and the acid OH
    assert atoms[1] == [2]
    assert bonds[1] == []  # a single atom has no internal bonds


def test_diff_highlights_show_added_atoms_only():
    ref = Chem.MolFromSmiles(ASPIRIN)
    mols = [
        Chem.MolFromSmiles(s)
        for s in (ASPIRIN, SALICYLIC, METHYL_SALICYLATE, "C1CCCCC1")
    ]
    atoms, _, not_computed = drawing.diff_highlights(mols, ref)
    assert atoms[0] == []  # identical to the reference
    assert atoms[1] == []  # a substructure of the reference adds nothing
    assert len(atoms[2]) == 1  # the ester methyl
    assert mols[2].GetAtomWithIdx(atoms[2][0]).GetSymbol() == "C"
    assert atoms[3] == []  # nothing in common: no highlight rather than everything
    assert not_computed == []


# FindMCS returns whatever it has when its 2 s run out, and the code took that for
# the maximum common substructure: dolutegravir against cabotegravir came back with
# a 2-atom "core" and 28 of 30 atoms highlighted as changed. About one analogue pair
# in eleven (Tanimoto >= 0.5 in the demo set) hits the limit.

DOLUTEGRAVIR = "C[C@@H]1CCO[C@@H]2N1C(=O)C3=C(C(=O)C(=CN3C2)C(=O)NCC4=C(C=C(C=C4)F)F)O"
CABOTEGRAVIR = "C[C@H]1CO[C@H]2N1C(=O)C3=C(C(=O)C(=CN3C2)C(=O)NCC4=C(C=C(C=C4)F)F)O"


class _CutShort:
    """What FindMCS returned for the pair above after 2 s."""

    canceled = True
    numAtoms = 2
    smartsString = "[#6]-[#6]"


def test_a_search_that_timed_out_is_not_taken_for_the_answer(monkeypatch):
    monkeypatch.setattr(drawing.rdFMCS, "FindMCS", lambda *a, **k: _CutShort())
    mol = Chem.MolFromSmiles(DOLUTEGRAVIR)
    atoms, bonds, not_computed = drawing.diff_highlights(
        [mol], Chem.MolFromSmiles(CABOTEGRAVIR)
    )
    # The ring that has one more atom than cabotegravir's, and its methyl: the 5
    # atoms the search finds when it is left to finish (30 s), not 28.
    assert sorted(mol.GetAtomWithIdx(a).GetSymbol() for a in atoms[0]) == [
        "C",
        "C",
        "C",
        "C",
        "O",
    ]
    assert bonds[0]
    assert not_computed == []


def test_the_fallback_agrees_with_a_search_that_finishes(monkeypatch):
    """Close analogues the 2 s search does finish: both give the same atoms."""
    pairs = [
        (ASPIRIN, METHYL_SALICYLATE),
        ("CC(C)Cc1ccc(cc1)C(C)C(=O)O", "CC(C)Cc1ccc(cc1)C(C)C(=O)OC"),
        ("c1ccc2ccccc2c1CCN", "c1ccc2ccccc2c1CCNC(C)=O"),
        ("O=C(O)c1ccccc1O", "O=C(O)c1ccc(Cl)cc1O"),
    ]
    finished = [
        drawing.diff_highlights([Chem.MolFromSmiles(b)], Chem.MolFromSmiles(a))[0]
        for a, b in pairs
    ]
    monkeypatch.setattr(drawing.rdFMCS, "FindMCS", lambda *a, **k: _CutShort())
    fallback = [
        drawing.diff_highlights([Chem.MolFromSmiles(b)], Chem.MolFromSmiles(a))[0]
        for a, b in pairs
    ]
    assert fallback == finished
    assert all(len(atoms[0]) >= 1 for atoms in finished)


def test_a_difference_that_cannot_be_computed_is_labelled(monkeypatch):
    monkeypatch.setattr(drawing.rdFMCS, "FindMCS", lambda *a, **k: _CutShort())
    monkeypatch.setattr(drawing.rdRascalMCES, "FindMCES", lambda *a, **k: [])
    seen = {}
    real = drawing.Draw.MolsToGridImage

    def spy(mols, **kwargs):
        seen.update(kwargs)
        return real(mols, **kwargs)

    monkeypatch.setattr(drawing.Draw, "MolsToGridImage", spy)
    png = drawing.draw_grid(
        [DOLUTEGRAVIR, DOLUTEGRAVIR],
        legends=["dolutegravir", ""],
        highlight="diff",
        reference=CABOTEGRAVIR,
    )
    assert png.startswith(b"\x89PNG")
    assert seen["legends"] == [
        "dolutegravir (difference not computed)",
        "difference not computed",
    ]
    assert seen["highlightAtomLists"] == [[], []]  # not the whole molecule


def test_atom_contributions_identical_molecule_is_all_non_negative():
    mol = Chem.MolFromSmiles(ASPIRIN)
    weights = drawing.atom_contributions(mol, mol)
    assert len(weights) == mol.GetNumAtoms()
    assert all(w >= 0 for w in weights)
    assert any(w > 0 for w in weights)


@pytest.mark.parametrize(
    ("kwargs", "error", "match"),
    [
        ({"smiles": []}, InvalidInput, "between 1 and 24"),
        ({"smiles": ["C"] * 25}, InvalidInput, "between 1 and 24"),
        ({"smiles": ["C"], "size": 100}, InvalidInput, "size"),
        ({"smiles": ["C"], "per_row": 9}, InvalidInput, "per_row"),
        ({"smiles": ["C", "CC"], "legends": ["one"]}, InvalidInput, "one legend"),
        ({"smiles": ["C"], "highlight": "diff"}, InvalidInput, "reference"),
        ({"smiles": ["C"], "highlight": "similarity_map"}, InvalidInput, "reference"),
        ({"smiles": ["C"], "highlight": "smarts"}, InvalidInput, "SMARTS"),
        (
            {"smiles": ["C"], "highlight": "smarts", "smarts": "[[["},
            InvalidInput,
            "Invalid SMARTS",
        ),
        ({"smiles": ["C"], "highlight": "sparkles"}, InvalidInput, "Unknown highlight"),
        ({"smiles": ["C", "nope("]}, MoleculeValidationError, "position 2"),
        (
            {"smiles": ["C"], "highlight": "diff", "reference": "nope("},
            MoleculeValidationError,
            "Invalid SMILES",
        ),
    ],
)
def test_validation(kwargs, error, match):
    with pytest.raises(error, match=match):
        drawing.draw_grid(**kwargs)


def test_oversized_molecules_and_patterns_are_refused():
    with pytest.raises(MoleculeValidationError, match="position 2 has 1,001 atoms"):
        drawing.draw_grid(["CCO", "C" * 1_001])
    with pytest.raises(MoleculeValidationError, match="position 1 is 10,001"):
        drawing.draw_grid(["C" * 10_001])
    with pytest.raises(MoleculeValidationError, match="SMARTS has 1,001 atoms"):
        drawing.draw_grid(["CCO"], highlight="smarts", smarts="[#6]" * 1_001)


def test_similarity_maps_are_bounded(monkeypatch):
    # Hardening audit C3: aligning a 606-atom chain to the core made its similarity
    # map take 41 s with the GIL held, freezing an HTTP server for every client.
    big = "c1ccccc1" + "CCO" * 150  # 456 atoms
    aligned: list[int] = []
    real = drawing.diff_highlights

    def spy(mols, reference):
        aligned.extend(m.GetNumAtoms() for m in mols)
        return real(mols, reference)

    monkeypatch.setattr(drawing, "diff_highlights", spy)
    png = drawing.draw_grid(
        ["c1ccccc1CCO", big], highlight="similarity_map", reference="c1ccccc1CO"
    )
    assert png.startswith(b"\x89PNG")
    assert aligned == [9]  # the long chain keeps RDKit's own layout
    with pytest.raises(InvalidInput, match="2,700 atoms in one request"):
        drawing.draw_grid(["C" * 900] * 3, highlight="similarity_map", reference="CCO")
