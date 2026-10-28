"""The version and registry name live in several files; they must agree."""

import json
import re
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


def test_icon_urls_point_at_a_file_in_the_repo():
    """The server and the registry entry show the plugin's icon from master; moving
    the file would break both silently."""
    from pubchem_finder.mcp_server import ICON_URL

    prefix = "https://raw.githubusercontent.com/MysterionRise/pubchem-finder/master/"
    urls = {ICON_URL, *(icon["src"] for icon in _json("server.json")["icons"])}
    assert urls == {ICON_URL}
    assert ICON_URL.startswith(prefix)
    assert (ROOT / ICON_URL.removeprefix(prefix)).is_file()


PLUGIN = ROOT / "plugins" / "pubchem-finder"


def _release_tuple(version: str) -> tuple[int, ...]:
    return tuple(int(part) for part in version.split("."))


def test_plugin_runs_one_released_version_from_its_lockfile():
    """Anthropic's plugin directory wants the server run from a lockfile with exact
    versions, registry sources and hashes. The plugin pins one released version in
    five places (scripts/plugin_release.py moves them together); it may lag the
    package between a release and that bump, but never run ahead of it."""
    pyproject = (PLUGIN / "pyproject.toml").read_text(encoding="utf-8")
    pins = re.findall(r'"pubchem-finder==([^"]+)"', pyproject)
    assert len(pins) == 1
    pinned = pins[0]

    lock = (PLUGIN / "uv.lock").read_text(encoding="utf-8")
    locked = re.search(
        r'^name = "pubchem-finder"\nversion = "([^"]+)"\nsource = \{ registry = "https://pypi.org/simple" \}',
        lock,
        re.M,
    )
    assert locked is not None, "pubchem-finder is not locked from PyPI"
    assert locked.group(1) == pinned
    assert 'hash = "sha256:' in lock
    assert (PLUGIN / "uv.lock").stat().st_size < 256 * 1024, (
        "over the directory's limit"
    )

    assert (
        _json("plugins/pubchem-finder/.claude-plugin/plugin.json")["version"] == pinned
    )
    assert _json(".claude-plugin/marketplace.json")["plugins"][0]["version"] == pinned
    assert _release_tuple(pinned) <= _release_tuple(__version__)


def test_plugin_command_runs_the_pinned_release():
    """The directory can pin and verify exactly `uvx <package>==<version>`: no other
    options, no environment."""
    pyproject = (PLUGIN / "pyproject.toml").read_text(encoding="utf-8")
    (pinned,) = re.findall(r'"pubchem-finder==([^"]+)"', pyproject)
    server = _json("plugins/pubchem-finder/.mcp.json")["mcpServers"]["pubchem-finder"]
    assert server == {"command": "uvx", "args": [f"pubchem-finder=={pinned}", "mcp"]}
