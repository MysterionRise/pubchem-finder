"""scripts/demo_media.py: collect results from an index, then render the animation."""

import importlib.util
import sys
from pathlib import Path

from PIL import Image

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "demo_media.py"


def load_script():
    spec = importlib.util.spec_from_file_location("demo_media", SCRIPT)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module  # dataclasses look their module up here
    spec.loader.exec_module(module)
    return module


def test_collect_then_render_gif(substructure_home, tmp_path):
    demo = load_script()
    data = demo.collect(
        substructure_home,
        "CC(=O)Oc1ccccc1C(=O)O",
        "aspirin",
        "[CX3](=O)[OX2H1]",
        "carboxylic acid",
    )
    assert data["similarity"]["hits"][0]["cid"] == 2244
    assert data["substructure"]["hits"]
    assert data["similarity_ms"] > 0
    assert data["substructure_ms"] > 0

    written = demo.render(
        data, tmp_path, example="OC(=O)c1ccccc1", fps=2, use_ffmpeg=False
    )
    with Image.open(written["gif"]) as gif:
        assert gif.size == (960, 540)
        assert gif.n_frames > 5  # identical consecutive frames are merged
        total_ms = 0
        for i in range(gif.n_frames):
            gif.seek(i)
            total_ms += gif.info["duration"]
    assert 12_000 < total_ms < 30_000


def test_single_component_skips_salts_and_mixtures():
    demo = load_script()
    hits = [
        {"cid": 1, "smiles": "CCO"},
        {"cid": 2, "smiles": "CC(=O)[O-].[Na+]"},
        {"cid": 3, "smiles": None},
        {"cid": 4, "smiles": "c1ccccc1"},
    ]
    assert [h["cid"] for h in demo.single_component(hits, 5)] == [1, 4]
    assert [h["cid"] for h in demo.single_component(hits, 1)] == [1]
