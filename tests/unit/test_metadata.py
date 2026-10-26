"""The version and registry name live in several files; they must agree."""

import json
from pathlib import Path

from pubchem_finder import __version__

ROOT = Path(__file__).resolve().parents[2]


def _json(relative: str) -> dict:
    return json.loads((ROOT / relative).read_text(encoding="utf-8"))


def test_versions_match_the_package():
    server = _json("server.json")
    found = {
        "server.json": server["version"],
        **{
            f"server.json packages[{i}]": package["version"]
            for i, package in enumerate(server["packages"])
        },
        "plugin.json": _json("plugins/pubchem-finder/.claude-plugin/plugin.json")[
            "version"
        ],
        "marketplace.json": _json(".claude-plugin/marketplace.json")["plugins"][0][
            "version"
        ],
    }
    cff = (ROOT / "CITATION.cff").read_text(encoding="utf-8").splitlines()
    found["CITATION.cff"] = next(
        line.split(":", 1)[1].strip() for line in cff if line.startswith("version:")
    )
    assert found == dict.fromkeys(found, __version__)


def test_registry_name_marker_is_in_the_readme():
    """The MCP registry proves PyPI ownership by finding `mcp-name: <name>` on its
    own line in the package description, which is the README."""
    name = _json("server.json")["name"]
    lines = (ROOT / "README.md").read_text(encoding="utf-8").splitlines()
    assert f"<!-- mcp-name: {name} -->" in lines
