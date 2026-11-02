"""What the index messages tell people to run, and what it costs.

The messages reach people where they use pubchem-finder: often an MCP client that
launched the server with uvx, with no ``pubchem-finder`` on their PATH and without
the server's environment. So commands run this exact release through uvx and name
the data directory unless it is the default.
"""

import re
from pathlib import Path

from . import __version__
from .config import default_home

LAUNCH_PYTHON = "3.14"
"""The Python every published launch command passes to uvx: FPSim2, numexpr and
ndindex publish no wheels for newer versions yet."""

BUILD_SAMPLE = "index build --limit 1 --with-substructure"
SAMPLE_COST = (
    "PubChem's first SDF chunk, about 445,000 compounds: a 323 MB download and "
    "about 5 minutes"
)
ALL_COST = (
    "all of PubChem, about 120 million compounds in 359 chunks: 111 GB to "
    "download, an overnight job"
)
_RELEASE = re.compile(r"\d+(\.\d+)*")


def command(args: str, home: Path | None = None) -> str:
    """A command line running ``pubchem-finder <args>`` from any terminal."""
    if _RELEASE.fullmatch(__version__) and __version__ != "0.0.0":
        program = f"uvx --python {LAUNCH_PYTHON} pubchem-finder=={__version__}"
    else:  # a development build: run what is installed
        program = "pubchem-finder"
    if home is not None and home != default_home():
        args += f" --home {_quote(str(home))}"
    return f"{program} {args}"


def build_index(home: Path | None = None) -> str:
    """How to build a first index, and what it costs."""
    return (
        f"Build one from a terminal with `{command(BUILD_SAMPLE, home)}` "
        f"({SAMPLE_COST}); `--all` instead of `--limit 1` indexes {ALL_COST}."
    )


def add_substructure(home: Path | None = None) -> str:
    """How to add (or refresh) substructure search on an existing index."""
    return (
        f"`{command('index build --with-substructure', home)}` (uses the compounds "
        "already indexed; no download)"
    )


def _quote(text: str) -> str:
    # Double quotes work in POSIX shells, PowerShell and cmd.exe alike.
    return f'"{text}"' if any(c.isspace() for c in text) else text
