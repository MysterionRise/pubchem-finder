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


LAUNCH_PYTHON = "3.14"


def test_published_launch_commands_pin_a_supported_python():
    """uvx runs uv's newest Python whatever the package supports, and FPSim2, numexpr
    and ndindex publish no wheels for 3.15: every command we publish pins 3.14."""
    import tomllib

    package = _json("server.json")["packages"][0]
    runtime = {(a["name"], a["value"]) for a in package["runtimeArguments"]}
    assert ("--python", LAUNCH_PYTHON) in runtime
    pinned = f'"args": ["--python", "{LAUNCH_PYTHON}", "pubchem-finder", "mcp"]'
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    clients = (ROOT / "docs" / "clients.md").read_text(encoding="utf-8")
    general = clients.split("\n## Clients")[
        0
    ]  # tested setups stay as their testers ran them
    for name, text in {"README.md": readme, "docs/clients.md": general}.items():
        assert pinned in text, name
        assert "uvx pubchem-finder" not in text, name
        assert '"args": ["pubchem-finder", "mcp"]' not in text, name
    example = (ROOT / "examples" / "mcp_client.py").read_text(encoding="utf-8")
    assert f"uvx --python {LAUNCH_PYTHON} pubchem-finder" in example
    project = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    assert (
        f"Programming Language :: Python :: {LAUNCH_PYTHON}"
        in (project["project"]["classifiers"])
    )


def test_plugin_release_moves_every_pin(tmp_path):
    """scripts/plugin_release.py finds the pinned argument wherever it sits, so a
    flag before it (such as --python) is never overwritten."""
    import importlib.util
    import shutil

    shutil.copytree(ROOT / "plugins", tmp_path / "plugins")
    (tmp_path / ".claude-plugin").mkdir()
    shutil.copy(
        ROOT / ".claude-plugin" / "marketplace.json", tmp_path / ".claude-plugin"
    )
    mcp_json = tmp_path / "plugins" / "pubchem-finder" / ".mcp.json"
    config = json.loads(mcp_json.read_text(encoding="utf-8"))
    server = config["mcpServers"]["pubchem-finder"]
    server["args"] = ["--python=3.14", "pubchem-finder==1.0.0", "mcp"]
    mcp_json.write_text(json.dumps(config), encoding="utf-8")
    calls: list[list[str]] = []
    spec = importlib.util.spec_from_file_location(
        "plugin_release", ROOT / "scripts" / "plugin_release.py"
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    module.main("9.9.9", root=tmp_path, run=lambda cmd, cwd, check: calls.append(cmd))

    plugin = tmp_path / "plugins" / "pubchem-finder"
    args = json.loads(mcp_json.read_text(encoding="utf-8"))["mcpServers"][
        "pubchem-finder"
    ]["args"]
    assert args == ["--python=3.14", "pubchem-finder==9.9.9", "mcp"]
    assert '"pubchem-finder==9.9.9"' in (plugin / "pyproject.toml").read_text()
    manifest = json.loads((plugin / ".claude-plugin" / "plugin.json").read_text())
    assert manifest["version"] == "9.9.9"
    market = json.loads((tmp_path / ".claude-plugin" / "marketplace.json").read_text())
    assert market["plugins"][0]["version"] == "9.9.9"
    plugin_dir = str(Path("plugins/pubchem-finder"))
    assert calls == [["uv", "lock", "--upgrade", "--project", plugin_dir]]


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
    """The directory pins and verifies `uvx <package>==<version>`; it blocks an `env`
    block and any variable but ${CLAUDE_PLUGIN_ROOT}. The one option keeps uv on a
    Python every dependency has wheels for (see test_published_launch_commands...)."""
    pyproject = (PLUGIN / "pyproject.toml").read_text(encoding="utf-8")
    (pinned,) = re.findall(r'"pubchem-finder==([^"]+)"', pyproject)
    server = _json("plugins/pubchem-finder/.mcp.json")["mcpServers"]["pubchem-finder"]
    assert server == {
        "command": "uvx",
        "args": [f"--python={LAUNCH_PYTHON}", f"pubchem-finder=={pinned}", "mcp"],
    }
