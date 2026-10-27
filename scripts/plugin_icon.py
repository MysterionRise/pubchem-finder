"""Render the Claude plugin icon: aspirin drawn in white on the README's blue.

Anthropic's plugin directory wants a square PNG of 512 to 2048 px under 2 MB at
``plugins/pubchem-finder/.claude-plugin/icon.png``. Run from the repository root:

    uv run python scripts/plugin_icon.py
"""

from __future__ import annotations

import io
import sys
from pathlib import Path

from PIL import Image
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

SIZE = 1024
BLUE = (0x1B / 255, 0x6F / 255, 0xA8 / 255, 1.0)
ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"
DEFAULT_OUT = Path("plugins/pubchem-finder/.claude-plugin/icon.png")


def render(out: Path = DEFAULT_OUT) -> Path:
    mol = Chem.MolFromSmiles(ASPIRIN)
    rdDepictor.SetPreferCoordGen(True)
    rdDepictor.Compute2DCoords(mol)

    drawer = rdMolDraw2D.MolDraw2DCairo(SIZE, SIZE)
    options = drawer.drawOptions()
    options.setBackgroundColour(BLUE)
    options.useBWAtomPalette()
    options.updateAtomPalette({-1: (1.0, 1.0, 1.0)})
    options.padding = 0.16  # type: ignore[assignment]  # stub says MolDrawOptions
    options.bondLineWidth = 7  # type: ignore[assignment]  # stub says MolDrawOptions
    options.scaleBondWidth = False  # type: ignore[assignment]  # stub says MolDrawOptions
    options.fixedFontSize = 110  # type: ignore[assignment]  # stub says MolDrawOptions
    options.additionalAtomLabelPadding = 0.12  # type: ignore[assignment]  # stub says MolDrawOptions
    options.multipleBondOffset = 0.18  # type: ignore[assignment]  # stub says MolDrawOptions
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()

    image = Image.open(io.BytesIO(drawer.GetDrawingText())).convert("RGB")
    out.parent.mkdir(parents=True, exist_ok=True)
    image.save(out, optimize=True)
    return out


if __name__ == "__main__":
    path = render(Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_OUT)
    print(f"{path}: {SIZE}x{SIZE} PNG, {path.stat().st_size / 1024:.0f} KB")
