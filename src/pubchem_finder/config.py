"""Where pubchem-finder keeps its data.

The data directory defaults to the platform's user data directory
(``~/.local/share/pubchem-finder`` on Linux, ``~/Library/Application Support/
pubchem-finder`` on macOS, ``%LOCALAPPDATA%\\pubchem-finder`` on Windows) and can
be overridden with ``PUBCHEM_FINDER_HOME`` or an explicit path.
"""

import os
from pathlib import Path

from platformdirs import user_data_dir

HOME_ENV = "PUBCHEM_FINDER_HOME"
INDEX_MODE_ENV = "PUBCHEM_FINDER_INDEX_MODE"
CACHE_ENV = "PUBCHEM_FINDER_CACHE"


def resolve_home(explicit: Path | str | None = None) -> Path:
    """Return the data directory: explicit path, then ``$PUBCHEM_FINDER_HOME``,
    then the platform default. The directory is not created."""
    if explicit is not None:
        return Path(explicit).expanduser().resolve()
    env = os.environ.get(HOME_ENV)
    if env:
        return Path(env).expanduser().resolve()
    return Path(user_data_dir("pubchem-finder", appauthor=False))


def cache_enabled() -> bool:
    """The PubChem response cache is on unless ``PUBCHEM_FINDER_CACHE`` is 0/false/off."""
    return os.environ.get(CACHE_ENV, "1").strip().lower() not in {
        "0",
        "false",
        "no",
        "off",
    }


def cache_path(home: Path | str | None = None) -> Path:
    """Location of the PubChem response cache inside the data directory."""
    return resolve_home(home) / "cache" / "pubchem-http.sqlite"
