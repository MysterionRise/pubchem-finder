"""Point the Claude plugin at a released pubchem-finder version and re-lock it.

Run after a release, once the version is on PyPI (the release workflow's
plugin-branch job does, and pushes the result for a PR)::

    uv run python scripts/plugin_release.py 2.1.0

It sets the pin in plugins/pubchem-finder/pyproject.toml and .mcp.json, the plugin and
marketplace versions, and re-locks plugins/pubchem-finder/uv.lock with the newest
dependencies, which is what `uvx pubchem-finder==<version>` installs.
tests/unit/test_metadata.py checks that the five agree.
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
from collections.abc import Callable
from pathlib import Path
from typing import Any

PLUGIN = Path("plugins/pubchem-finder")


def _edit_json(path: Path, key_path: list[str | int], value: str) -> None:
    data = json.loads(path.read_text(encoding="utf-8"))
    target = data
    for key in key_path[:-1]:
        target = target[key]
    target[key_path[-1]] = value
    path.write_text(
        json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )


def main(
    version: str,
    root: Path = Path("."),
    run: Callable[..., Any] = subprocess.run,
) -> None:
    """Point the plugin under ``root`` (the repository) at ``version`` and re-lock."""
    if not re.fullmatch(r"\d+\.\d+\.\d+", version):
        sys.exit(f"not a release version: {version}")
    plugin = root / PLUGIN
    pyproject = plugin / "pyproject.toml"
    text, count = re.subn(
        r'"pubchem-finder==[^"]+"',
        f'"pubchem-finder=={version}"',
        pyproject.read_text(encoding="utf-8"),
    )
    if count != 1:
        sys.exit(f"expected one pubchem-finder pin in {pyproject}, found {count}")
    mcp_json = plugin / ".mcp.json"
    args = json.loads(mcp_json.read_text(encoding="utf-8"))["mcpServers"][
        "pubchem-finder"
    ]["args"]
    pins = [i for i, arg in enumerate(args) if arg.startswith("pubchem-finder==")]
    if len(pins) != 1:
        sys.exit(f"expected one pubchem-finder== argument in {mcp_json}, found {args}")
    pyproject.write_text(text, encoding="utf-8")
    _edit_json(plugin / ".claude-plugin" / "plugin.json", ["version"], version)
    _edit_json(
        mcp_json,
        ["mcpServers", "pubchem-finder", "args", pins[0]],
        f"pubchem-finder=={version}",
    )
    _edit_json(
        root / ".claude-plugin" / "marketplace.json", ["plugins", 0, "version"], version
    )
    run(["uv", "lock", "--upgrade", "--project", str(PLUGIN)], cwd=root, check=True)
    print(f"plugin now runs pubchem-finder {version}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
