"""Render the search walkthrough animation (README hero GIF, MP4 for posts).

Two steps, so the animation shows real results and real timings:

    # 1. run the queries against your index; saves hits and warm query times
    uv run python scripts/demo_media.py collect --home DATA_DIR \\
        --query "CC(=O)Oc1ccccc1C(=O)O" --query-name aspirin \\
        --smarts "[NX3][SX4](=O)(=O)c" --smarts-name "aryl sulfonamide" -o demo-data.json
    # 2. draw it
    uv run python scripts/demo_media.py render demo-data.json -o docs/media

render writes search-demo.mp4 (1280x720, 30 fps) and search-demo.gif (960x540). ffmpeg
makes both; without it only a GIF is written, with Pillow.
"""

from __future__ import annotations

import argparse
import io
import json
import shutil
import statistics
import subprocess
import time
from collections.abc import Callable, Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from PIL import Image, ImageDraw, ImageFont
from rdkit import Chem, rdBase
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

from pubchem_finder.drawing import MIN_FONT_SIZE, diff_highlights, smarts_highlights

rdBase.DisableLog("rdApp.*")

W, H = 1280, 720
GIF_WIDTH = 960
ACCENT = (16, 110, 160)
INK = (20, 24, 31)
MUTED = (100, 110, 125)
DIFF_COLOR = (1.0, 0.72, 0.42)  # warm: what changed relative to the query
MATCH_COLOR = (0.55, 0.78, 1.0)  # cool: the substructure match
REPO = "github.com/MysterionRise/pubchem-finder"

FONT_DIRS = [Path("/usr/share/fonts/opentype/inter"), Path("/Library/Fonts")]
MONO = Path("/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf")


def font(weight: str, size: int) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    """Inter if installed, else Pillow's built-in font (same layout, plainer look)."""
    for directory in FONT_DIRS:
        path = directory / f"Inter-{weight}.otf"
        if path.exists():
            return ImageFont.truetype(str(path), size)
    if weight == "Mono" and MONO.exists():
        return ImageFont.truetype(str(MONO), size)
    return ImageFont.load_default(size)


def fit(text: str, face: Any, width: int) -> str:
    """Truncate ``text`` with an ellipsis so it fits in ``width`` pixels."""
    if face.getlength(text) <= width:
        return text
    while text and face.getlength(text + "…") > width:
        text = text[:-1]
    return text.rstrip(" -,;") + "…"


def draw_mol(
    mol: Chem.Mol,
    size: tuple[int, int],
    atoms: list[int] | None = None,
    bonds: list[int] | None = None,
    color: tuple[float, float, float] = DIFF_COLOR,
) -> Image.Image:
    drawer = rdMolDraw2D.MolDraw2DCairo(*size)
    options = drawer.drawOptions()
    options.minFontSize = MIN_FONT_SIZE  # type: ignore[assignment]  # stub says MolDrawOptions
    options.clearBackground = True  # type: ignore[assignment]  # stub says MolDrawOptions
    atoms, bonds = atoms or [], bonds or []
    rdMolDraw2D.PrepareAndDrawMolecule(
        drawer,
        mol,
        highlightAtoms=atoms,
        highlightBonds=bonds,
        highlightAtomColors=dict.fromkeys(atoms, color),
        highlightBondColors=dict.fromkeys(bonds, color),
    )
    drawer.FinishDrawing()
    return Image.open(io.BytesIO(drawer.GetDrawingText())).convert("RGB")


@dataclass
class Tile:
    image: Image.Image
    title: str
    subtitle: str
    detail: str


TILE_W, TILE_H = 200, 290
GRID_X, GRID_Y = 410, 100


def render_tile(tile: Tile) -> Image.Image:
    card = Image.new("RGB", (TILE_W, TILE_H), "white")
    card.paste(tile.image, (2, 2))
    d = ImageDraw.Draw(card)
    d.rounded_rectangle(
        [0, 0, TILE_W - 1, TILE_H - 1], radius=10, outline=(225, 229, 235)
    )
    bold, regular, small = (
        font("SemiBold", 17),
        font("Regular", 15),
        font("Regular", 13),
    )
    d.text((12, 212), tile.title, font=bold, fill=INK)
    d.text((12, 236), tile.subtitle, font=regular, fill=ACCENT)
    d.text((12, 260), fit(tile.detail, small, TILE_W - 24), font=small, fill=MUTED)
    return card


def single_component(hits: list[dict[str, Any]], n: int) -> list[dict[str, Any]]:
    """The first ``n`` hits that are one molecule (no salts or mixtures)."""
    return [h for h in hits if h.get("smiles") and "." not in h["smiles"]][:n]


def similarity_tiles(search: dict[str, Any], n: int = 8) -> tuple[Chem.Mol, list[Tile]]:
    query = Chem.MolFromSmiles(search["query_smiles"])
    hits = single_component(search["hits"], n)
    mols = [Chem.MolFromSmiles(h["smiles"]) for h in hits]
    atoms, bonds = diff_highlights(mols, query)  # also aligns depictions to the query
    tiles = [
        Tile(
            draw_mol(mol, (TILE_W - 4, 205), a, b, DIFF_COLOR),
            f"CID {h['cid']}",
            f"Tanimoto {h['similarity']:.3f}",
            h.get("iupac_name") or "",
        )
        for h, mol, a, b in zip(hits, mols, atoms, bonds, strict=True)
    ]
    return query, tiles


def substructure_tiles(sub: dict[str, Any], n: int = 8) -> list[Tile]:
    hits = single_component(sub["hits"], n)
    mols = [Chem.MolFromSmiles(h["smiles"]) for h in hits]
    for mol in mols:
        rdDepictor.Compute2DCoords(mol)
    atoms, bonds = smarts_highlights(mols, sub["query"])
    return [
        Tile(
            draw_mol(mol, (TILE_W - 4, 205), a, b, MATCH_COLOR),
            f"CID {h['cid']}",
            h.get("molecular_formula") or "",
            h.get("iupac_name") or "",
        )
        for h, mol, a, b in zip(hits, mols, atoms, bonds, strict=True)
    ]


def base(subtitle: str) -> Image.Image:
    frame = Image.new("RGB", (W, H), "white")
    d = ImageDraw.Draw(frame)
    d.rectangle([0, 0, W, 6], fill=ACCENT)
    d.text((40, 24), "pubchem-finder", font=font("Bold", 30), fill=INK)
    d.text((W - 40, 34), subtitle, font=font("Regular", 20), fill=MUTED, anchor="ra")
    return frame


def grid_position(i: int) -> tuple[int, int]:
    return GRID_X + (i % 4) * (TILE_W + 10), GRID_Y + (i // 4) * (TILE_H + 12)


def ease(t: float) -> float:
    t = min(max(t, 0.0), 1.0)
    return 1 - (1 - t) ** 3


def paste_faded(
    frame: Image.Image, card: Image.Image, xy: tuple[int, int], alpha: float
) -> None:
    if alpha <= 0:
        return
    x, y = xy
    y += int(18 * (1 - ease(alpha)))  # slide up a little while fading in
    region = frame.crop((x, y, x + card.width, y + card.height))
    frame.paste(Image.blend(region, card, ease(alpha)), (x, y))


Scene = Callable[[float], Image.Image]


def similarity_scene(
    query: Chem.Mol,
    query_name: str,
    tiles: list[Tile],
    search: dict[str, Any],
    ms: float | None,
) -> tuple[float, Scene]:
    cards = [render_tile(t) for t in tiles]
    query_img = draw_mol(query, (360, 300))
    n_indexed = f"{search['n_indexed']:,}"
    n = len(search["hits"])
    capped = "top " if n == search.get("top_k") else ""
    found = f"{capped}{n} hits ≥ {search['threshold']} Tanimoto"
    timing = f" in {ms:g} ms" if ms is not None else ""
    start, gap, duration = 1.6, 0.6, 1.6 + 0.6 * len(cards) + 1.6

    def scene(t: float) -> Image.Image:
        frame = base("similarity search over a local PubChem index")
        d = ImageDraw.Draw(frame)
        d.text((40, 96), "QUERY", font=font("SemiBold", 16), fill=MUTED)
        frame.paste(query_img, (30, 120))
        d.text((40, 424), query_name, font=font("Bold", 28), fill=INK)
        d.text((40, 462), search["query_smiles"], font=font("Mono", 15), fill=MUTED)
        sweep = ease((t - 0.3) / 0.6)
        d.rounded_rectangle([40, 510, 380, 518], radius=4, fill=(232, 236, 241))
        if sweep > 0:
            d.rounded_rectangle(
                [40, 510, 40 + int(340 * sweep), 518], radius=4, fill=ACCENT
            )
        d.text(
            (40, 530),
            f"Searching {n_indexed} PubChem compounds",
            font=font("Regular", 18),
            fill=INK,
        )
        if t >= 1.0:
            d.text(
                (40, 560), f"✓ {found}{timing}", font=font("SemiBold", 18), fill=ACCENT
            )
            d.text(
                (40, 594),
                "Highlighted: atoms that differ from the query.",
                font=font("Regular", 15),
                fill=MUTED,
            )
            d.text(
                (40, 616),
                "Showing single-component hits.",
                font=font("Regular", 15),
                fill=MUTED,
            )
        for i, card in enumerate(cards):
            paste_faded(frame, card, grid_position(i), (t - start - i * gap) / 0.35)
        return frame

    return duration, scene


def substructure_scene(
    tiles: list[Tile], sub: dict[str, Any], name: str, example: str, ms: float | None
) -> tuple[float, Scene]:
    cards = [render_tile(t) for t in tiles]
    example_mol = Chem.MolFromSmiles(example)
    atoms, bonds = smarts_highlights([example_mol], sub["query"])
    example_img = draw_mol(example_mol, (360, 250), atoms[0], bonds[0], MATCH_COLOR)
    timing = f" in {ms:g} ms" if ms is not None else ""
    start, gap, duration = 1.4, 0.55, 1.4 + 0.55 * len(cards) + 1.8

    def scene(t: float) -> Image.Image:
        frame = base("substructure search with SMILES or SMARTS")
        d = ImageDraw.Draw(frame)
        d.text((40, 96), "SUBSTRUCTURE", font=font("SemiBold", 16), fill=MUTED)
        d.text((40, 122), name, font=font("Bold", 28), fill=INK)
        d.text((40, 162), sub["query"], font=font("Mono", 18), fill=ACCENT)
        frame.paste(example_img, (30, 196))
        lines = [
            f"{sub['n_indexed']:,} compounds",
            f"→ {sub['n_candidates']:,} pass the fingerprint screen",
            f"→ checked atom by atom{timing}",
        ]
        for i, line in enumerate(lines):
            alpha = ease((t - 0.2 - 0.3 * i) / 0.3)
            if alpha > 0:
                shade = tuple(int(255 - (255 - c) * alpha) for c in INK)
                d.text((40, 470 + 30 * i), line, font=font("Regular", 19), fill=shade)
        if t >= 1.4:
            d.text(
                (40, 580),
                "Highlighted: the matched group.",
                font=font("Regular", 15),
                fill=MUTED,
            )
        for i, card in enumerate(cards):
            paste_faded(frame, card, grid_position(i), (t - start - i * gap) / 0.35)
        return frame

    return duration, scene


def end_scene() -> tuple[float, Scene]:
    def scene(t: float) -> Image.Image:
        frame = Image.new("RGB", (W, H), "white")
        d = ImageDraw.Draw(frame)
        d.rectangle([0, 0, W, 6], fill=ACCENT)
        d.text(
            (W / 2, 230), "pubchem-finder", font=font("Bold", 72), fill=INK, anchor="mm"
        )
        d.text(
            (W / 2, 310),
            "Search PubChem by structure, from your terminal or from Claude",
            font=font("Regular", 26),
            fill=MUTED,
            anchor="mm",
        )
        d.rounded_rectangle(
            [W / 2 - 250, 370, W / 2 + 250, 430], radius=10, fill=(24, 28, 36)
        )
        d.text(
            (W / 2, 400),
            "uvx pubchem-finder mcp",
            font=font("Mono", 26),
            fill="white",
            anchor="mm",
        )
        d.text(
            (W / 2, 480),
            "CLI  ·  Python library  ·  MCP server  ·  local similarity & substructure index",
            font=font("Regular", 20),
            fill=INK,
            anchor="mm",
        )
        d.text((W / 2, 530), REPO, font=font("SemiBold", 20), fill=ACCENT, anchor="mm")
        return frame

    return 2.5, scene


def crossfade(a: Scene, a_end: float, b: Scene, length: float = 0.4) -> Scene:
    def scene(t: float) -> Image.Image:
        return Image.blend(a(a_end), b(0.0), ease(t / length))

    return scene


def timeline(scenes: list[tuple[float, Scene]], fps: int) -> Iterator[Image.Image]:
    """Frames at ``fps``, with a short crossfade between scenes and back to the start."""
    fade = 0.4
    for i, (duration, scene) in enumerate(scenes):
        for k in range(round(duration * fps)):
            yield scene(k / fps)
        _, next_scene = scenes[(i + 1) % len(scenes)]
        blend = crossfade(scene, duration, next_scene, fade)
        for k in range(round(fade * fps)):
            yield blend(k / fps)


def render(
    data: dict[str, Any],
    out_dir: Path,
    example: str = "NS(=O)(=O)c1ccccc1",
    fps: int = 30,
    use_ffmpeg: bool = True,
) -> dict[str, Path]:
    """Render the animation from ``collect`` output; returns the paths written."""
    out_dir.mkdir(parents=True, exist_ok=True)
    search, substructure = data["similarity"], data["substructure"]
    query, sim_tiles = similarity_tiles(search)
    scenes = [
        similarity_scene(
            query, data["query_name"], sim_tiles, search, data.get("similarity_ms")
        ),
        substructure_scene(
            substructure_tiles(substructure),
            substructure,
            data["smarts_name"],
            example,
            data.get("substructure_ms"),
        ),
        end_scene(),
    ]
    ffmpeg = shutil.which("ffmpeg") if use_ffmpeg else None
    written: dict[str, Path] = {}
    if ffmpeg:
        mp4 = out_dir / "search-demo.mp4"
        encode = subprocess.Popen(
            [ffmpeg, "-y", "-loglevel", "error", "-f", "rawvideo", "-pix_fmt", "rgb24",
             "-s", f"{W}x{H}", "-r", str(fps), "-i", "-",
             "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "20", "-preset", "slow",
             "-movflags", "+faststart", str(mp4)],
            stdin=subprocess.PIPE,
        )  # fmt: skip
        assert encode.stdin is not None
        for frame in timeline(scenes, fps):
            encode.stdin.write(frame.tobytes())
        encode.stdin.close()
        if encode.wait() != 0:
            raise RuntimeError("ffmpeg failed to encode the MP4")
        gif = out_dir / "search-demo.gif"
        palette = (
            f"fps=12,scale={GIF_WIDTH}:-1:flags=lanczos,split[a][b];"
            "[a]palettegen=max_colors=128:stats_mode=diff[p];"
            "[b][p]paletteuse=dither=bayer:bayer_scale=5:diff_mode=rectangle"
        )
        subprocess.run(
            [ffmpeg, "-y", "-loglevel", "error", "-i", str(mp4), "-vf", palette, "-loop", "0", str(gif)],
            check=True,
        )  # fmt: skip
        written.update(mp4=mp4, gif=gif)
    else:
        gif = out_dir / "search-demo.gif"
        size = (GIF_WIDTH, GIF_WIDTH * H // W)
        frames = [
            f.resize(size, Image.Resampling.LANCZOS) for f in timeline(scenes, fps)
        ]
        frames[0].save(
            gif, save_all=True, append_images=frames[1:], duration=round(1000 / fps), loop=0,
            optimize=True,
        )  # fmt: skip
        written["gif"] = gif
    return written


def median_ms(run: Callable[[], object], repeats: int = 5) -> float:
    times = []
    for _ in range(repeats):
        start = time.perf_counter()
        run()
        times.append((time.perf_counter() - start) * 1000)
    return round(statistics.median(times), 1)


def collect(
    home: Path,
    query: str,
    query_name: str,
    smarts: str,
    smarts_name: str,
    threshold: float = 0.6,
    top: int = 40,
    max_hits: int = 24,
) -> dict[str, Any]:
    """Run both searches on the index in ``home``, with warm median timings
    measured through the same service the CLI and MCP server use."""
    from pubchem_finder import service

    def similar() -> Any:
        return service.similarity_search(
            query, threshold=threshold, top_k=top, home=home
        )

    def containing() -> Any:
        return service.substructure_search(
            smarts, query_type="smarts", max_hits=max_hits, home=home
        )

    similarity, substructure = similar(), containing()  # also loads the index
    return {
        "query_name": query_name,
        "smarts_name": smarts_name,
        "similarity": {**similarity.model_dump(mode="json"), "top_k": top},
        "substructure": substructure.model_dump(mode="json"),
        "similarity_ms": median_ms(similar),
        "substructure_ms": median_ms(containing),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    commands = parser.add_subparsers(dest="command", required=True)
    gather = commands.add_parser("collect", help="run the searches on an index")
    gather.add_argument("--home", type=Path, required=True, help="data directory")
    gather.add_argument("--query", required=True, help="similarity query SMILES")
    gather.add_argument("--query-name", required=True, help="caption for the query")
    gather.add_argument("--smarts", required=True, help="substructure query")
    gather.add_argument("--smarts-name", required=True, help="caption for it")
    gather.add_argument("-o", "--out", type=Path, default=Path("demo-data.json"))
    draw = commands.add_parser("render", help="draw the animation")
    draw.add_argument("data", type=Path, help="output of collect")
    draw.add_argument("-o", "--out", type=Path, default=Path("docs/media"))
    draw.add_argument(
        "--example",
        default="NS(=O)(=O)c1ccccc1",
        help="molecule that shows the pattern",
    )
    draw.add_argument(
        "--no-ffmpeg", action="store_true", help="write a Pillow GIF only"
    )
    args = parser.parse_args()

    if args.command == "collect":
        data = collect(
            args.home, args.query, args.query_name, args.smarts, args.smarts_name
        )
        args.out.write_text(json.dumps(data), encoding="utf-8")
        print(
            f"wrote {args.out}: {len(data['similarity']['hits'])} similar, "
            f"{len(data['substructure']['hits'])} containing; "
            f"{data['similarity_ms']} ms / {data['substructure_ms']} ms"
        )
        return
    written = render(
        json.loads(args.data.read_text(encoding="utf-8")),
        args.out,
        example=args.example,
        use_ffmpeg=not args.no_ffmpeg,
    )
    for path in written.values():
        print(f"wrote {path} ({path.stat().st_size / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()
