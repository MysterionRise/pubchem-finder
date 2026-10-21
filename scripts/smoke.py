"""Smoke-test an installed pubchem-finder: the CLI, then every MCP tool over stdio.

    python scripts/smoke.py                          # python -m pubchem_finder, offline checks
    python scripts/smoke.py --command .venv/bin/pubchem-finder --home DIR --live

Offline checks always run. ``--home`` adds the local-index tools (substructure
search too, if the index has it); ``--live`` adds the tools that call PubChem.
Prints a Markdown table (also appended to ``$GITHUB_STEP_SUMMARY`` when set) and
exits non-zero if any check fails. Needs the ``mcp`` package, which
pubchem-finder depends on, so run it with the interpreter it is installed in.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
import time
from collections.abc import Awaitable, Callable
from pathlib import Path
from typing import Any

import anyio
from mcp.client.session import ClientSession
from mcp.client.stdio import StdioServerParameters, stdio_client

N_TOOLS = 17
ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"
INDOLE = "c1ccc2[nH]ccc2c1"


class Smoke:
    def __init__(self) -> None:
        self.rows: list[tuple[str, bool, float, str]] = []

    def record(self, name: str, ok: bool, seconds: float, detail: str) -> None:
        self.rows.append((name, ok, seconds, detail))
        print(
            f"{'ok  ' if ok else 'FAIL'} {seconds:7.2f}s  {name}: {detail}", flush=True
        )

    def table(self) -> str:
        lines = ["| Check | Result | Seconds | Detail |", "|---|---|---|---|"]
        for name, ok, seconds, detail in self.rows:
            detail = detail.replace("|", "\\|")
            lines.append(
                f"| `{name}` | {'ok' if ok else '**FAIL**'} | {seconds:.2f} | {detail} |"
            )
        return "\n".join(lines)

    @property
    def failed(self) -> int:
        return sum(not ok for _, ok, _, _ in self.rows)


def run_cli(smoke: Smoke, command: list[str], env: dict[str, str]) -> None:
    def check(name: str, args: list[str], verify: Callable[[str], str]) -> None:
        start = time.perf_counter()
        try:
            out = subprocess.run(
                [*command, *args],
                env=env,
                capture_output=True,
                text=True,
                encoding="utf-8",
                timeout=300,
                check=True,
            ).stdout
            smoke.record(name, True, time.perf_counter() - start, verify(out))
        except Exception as exc:  # any failure is a failed check
            detail = getattr(exc, "stderr", None) or str(exc)
            smoke.record(
                name, False, time.perf_counter() - start, detail.strip()[-300:]
            )

    def props(out: str) -> str:
        data = json.loads(out)
        assert data["canonical_smiles"] == ASPIRIN, data
        return f"MW {data['molecular_weight']:.2f}"

    png = Path(tempfile.mkdtemp()) / "grid.png"

    def drawn(_: str) -> str:
        data = png.read_bytes()
        assert data.startswith(b"\x89PNG"), "not a PNG"
        return f"{len(data):,} byte PNG"

    check("cli --version", ["--version"], lambda out: out.strip())
    check("cli props", ["props", ASPIRIN, "--json"], props)
    check("cli draw", ["draw", ASPIRIN, "O=C(O)c1ccccc1O", "-o", str(png)], drawn)


Verify = Callable[[Any], str]


async def run_mcp(
    smoke: Smoke, command: list[str], env: dict[str, str], home: bool, live: bool
) -> None:
    params = StdioServerParameters(
        command=command[0], args=[*command[1:], "mcp"], env=env
    )
    async with (
        stdio_client(params) as (read, write),
        ClientSession(read, write) as session,
    ):

        async def check(
            name: str, call: Callable[[], Awaitable[Any]], verify: Verify
        ) -> Any:
            start = time.perf_counter()
            try:
                result = await call()
                if getattr(result, "is_error", False):
                    raise RuntimeError(result.content[0].text)
                smoke.record(name, True, time.perf_counter() - start, verify(result))
                return result
            except Exception as exc:
                smoke.record(name, False, time.perf_counter() - start, str(exc)[-300:])
                return None

        def tool(name: str, args: dict[str, Any]) -> Callable[[], Awaitable[Any]]:
            return lambda: session.call_tool(name, args)

        def data(result: Any) -> dict[str, Any]:
            return dict(result.structured_content)

        def expect(cond: bool, detail: str) -> str:
            if not cond:
                raise AssertionError(detail)
            return detail

        await check(
            "mcp initialize",
            session.initialize,
            lambda r: f"{r.server_info.name} {r.server_info.version}",
        )
        await check(
            "mcp list_tools",
            session.list_tools,
            lambda r: expect(
                len(r.tools) == N_TOOLS, f"{len(r.tools)} tools (expected {N_TOOLS})"
            ),
        )
        await check(
            "compound_properties",
            tool("compound_properties", {"smiles": ASPIRIN}),
            lambda r: expect(
                data(r)["canonical_smiles"] == ASPIRIN,
                f"MW {data(r)['molecular_weight']:.2f}",
            ),
        )
        await check(
            "standardize",
            tool("standardize", {"smiles": "CC(=O)[O-].[Na+]"}),
            lambda r: expect(data(r)["smiles"] == "CC(=O)O", data(r)["smiles"]),
        )
        await check(
            "compare",
            tool("compare", {"smiles_a": ASPIRIN, "smiles_b": "O=C(O)c1ccccc1O"}),
            lambda r: f"ECFP4 {data(r)['ecfp4_tanimoto']:.3f}",
        )
        await check(
            "medchem_profile",
            tool("medchem_profile", {"smiles": ASPIRIN}),
            lambda r: f"QED {data(r)['qed']:.2f}, SA {data(r)['sa_score']:.2f}",
        )
        await check(
            "draw_molecules",
            tool(
                "draw_molecules",
                {
                    "smiles": [ASPIRIN, "O=C(O)c1ccccc1O"],
                    "highlight": "diff",
                    "reference_smiles": ASPIRIN,
                },
            ),
            lambda r: expect(
                r.content[0].mime_type == "image/png", r.content[0].mime_type
            ),
        )
        status = await check(
            "index_status",
            tool("index_status", {}),
            lambda r: f"{data(r)['state']}, {data(r)['n_compounds']:,} compounds",
        )

        if home and status is not None:
            await check(
                "similarity_search",
                tool(
                    "similarity_search",
                    {"smiles": ASPIRIN, "threshold": 0.5, "top_k": 10},
                ),
                lambda r: expect(
                    data(r)["hits"][0]["cid"] == 2244,
                    f"{len(data(r)['hits'])} hits, first CID {data(r)['hits'][0]['cid']}",
                ),
            )
            if data(status).get("substructure"):
                for label, query, kind in (
                    ("indole", INDOLE, "smiles"),
                    ("aryl sulfonamide", "[NX3][SX4](=O)(=O)c", "smarts"),
                ):
                    await check(
                        f"substructure_search ({label})",
                        tool(
                            "substructure_search",
                            {"query": query, "query_type": kind, "max_hits": 20},
                        ),
                        lambda r: expect(
                            bool(data(r)["hits"]),
                            f"{len(data(r)['hits'])} hits of {data(r)['n_candidates']:,} candidates",
                        ),
                    )

        if live:
            await check(
                "pubchem_lookup",
                tool("pubchem_lookup", {"query": "caffeine"}),
                lambda r: expect(
                    data(r)["compounds"][0]["cid"] == 2519,
                    f"CID {data(r)['compounds'][0]['cid']}",
                ),
            )
            await check(
                "pubchem_similarity_search",
                tool(
                    "pubchem_similarity_search", {"smiles": ASPIRIN, "max_records": 10}
                ),
                lambda r: expect(bool(data(r)["hits"]), f"{len(data(r)['hits'])} hits"),
            )
            await check(
                "pubchem_substructure_search",
                tool(
                    "pubchem_substructure_search", {"smiles": INDOLE, "max_records": 10}
                ),
                lambda r: expect(bool(data(r)["hits"]), f"{len(data(r)['hits'])} hits"),
            )
            await check(
                "pubchem_exact_match",
                tool("pubchem_exact_match", {"smiles": ASPIRIN}),
                lambda r: expect(
                    any(h["cid"] == 2244 for h in data(r)["hits"]), "CID 2244"
                ),
            )
            await check(
                "pubchem_synonyms",
                tool("pubchem_synonyms", {"cid": 2244}),
                lambda r: expect(
                    "50-78-2" in data(r)["cas_numbers"],
                    f"CAS {', '.join(data(r)['cas_numbers'])}",
                ),
            )
            await check(
                "pubchem_xrefs",
                tool("pubchem_xrefs", {"cid": 2244}),
                lambda r: ", ".join(f"{k} {v:,}" for k, v in data(r)["counts"].items()),
            )
            await check(
                "pubchem_description",
                tool("pubchem_description", {"cid": 2244}),
                lambda r: expect(
                    bool(data(r)["descriptions"]),
                    f"{len(data(r)['descriptions'])} descriptions",
                ),
            )
            await check(
                "pubchem_annotations",
                tool(
                    "pubchem_annotations",
                    {"cid": 2244, "heading": "GHS Classification"},
                ),
                lambda r: expect(
                    bool(data(r)["annotations"]),
                    f"{len(data(r)['annotations'])} annotations",
                ),
            )
            await check(
                "pubchem_bioactivity",
                tool("pubchem_bioactivity", {"cid": 2244}),
                lambda r: expect(
                    data(r)["n_assays"] > 0, f"{data(r)['n_assays']:,} assays"
                ),
            )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument(
        "--command",
        help="pubchem-finder executable (default: python -m pubchem_finder)",
    )
    parser.add_argument("--home", type=Path, help="data directory with a built index")
    parser.add_argument(
        "--live", action="store_true", help="also call the PubChem tools"
    )
    args = parser.parse_args()

    command = (
        [args.command] if args.command else [sys.executable, "-m", "pubchem_finder"]
    )
    env = {**os.environ, "PUBCHEM_FINDER_CACHE": "0", "PYTHONIOENCODING": "utf-8"}
    if args.home:
        env["PUBCHEM_FINDER_HOME"] = str(args.home)
    else:
        env["PUBCHEM_FINDER_HOME"] = tempfile.mkdtemp(prefix="pubchem-finder-smoke-")

    smoke = Smoke()
    run_cli(smoke, command, env)
    anyio.run(run_mcp, smoke, command, env, args.home is not None, args.live)

    table = smoke.table()
    print(
        f"\n{table}\n\n{len(smoke.rows) - smoke.failed}/{len(smoke.rows)} checks passed"
    )
    if summary := os.environ.get("GITHUB_STEP_SUMMARY"):
        with open(summary, "a", encoding="utf-8") as fh:
            fh.write(f"### pubchem-finder smoke test\n\n{table}\n\n")
    return 1 if smoke.failed else 0


if __name__ == "__main__":
    sys.exit(main())
