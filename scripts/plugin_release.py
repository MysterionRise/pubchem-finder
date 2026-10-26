"""Point the Claude plugin at a released pubchem-finder version and re-lock it.

Run after a release, once the version is on PyPI::

    uv run python scripts/plugin_release.py 2.1.0

It sets the pin in plugins/pubchem-finder/pyproject.toml, the plugin and marketplace
versions, and refreshes plugins/pubchem-finder/uv.lock. tests/unit/test_metadata.py
checks that the four agree.
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

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


def main(version: str) -> None:
    if not re.fullmatch(r"\d+\.\d+\.\d+", version):
        sys.exit(f"not a release version: {version}")
    pyproject = PLUGIN / "pyproject.toml"
    text, count = re.subn(
        r'"pubchem-finder==[^"]+"',
        f'"pubchem-finder=={version}"',
        pyproject.read_text(encoding="utf-8"),
    )
    if count != 1:
        sys.exit(f"expected one pubchem-finder pin in {pyproject}, found {count}")
    pyproject.write_text(text, encoding="utf-8")
    _edit_json(PLUGIN / ".claude-plugin" / "plugin.json", ["version"], version)
    _edit_json(
        Path(".claude-plugin/marketplace.json"), ["plugins", 0, "version"], version
    )
    subprocess.run(["uv", "lock", "--project", str(PLUGIN)], check=True)
    print(f"plugin now runs pubchem-finder {version}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
