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
CACHE_MAX_MB_ENV = "PUBCHEM_FINDER_CACHE_MAX_MB"
DEMO_ENV = "PUBCHEM_FINDER_DEMO"


def resolve_home(explicit: Path | str | None = None) -> Path:
    """Return the data directory: explicit path, then ``$PUBCHEM_FINDER_HOME``,
    then the platform default. The directory is not created."""
    if explicit is not None:
        return Path(explicit).expanduser().resolve()
    env = os.environ.get(HOME_ENV)
    if env:
        return Path(env).expanduser().resolve()
    return default_home()


def default_home() -> Path:
    """The platform's user data directory for pubchem-finder."""
    return Path(user_data_dir("pubchem-finder", appauthor=False))


def cache_enabled() -> bool:
    """The PubChem response cache is on unless ``PUBCHEM_FINDER_CACHE`` is 0/false/off."""
    return _flag(CACHE_ENV)


def demo_enabled() -> bool:
    """Local searches fall back to the built-in demo index when no PubChem index is
    built, unless ``PUBCHEM_FINDER_DEMO`` is 0/false/off."""
    return _flag(DEMO_ENV)


def _flag(name: str) -> bool:
    return os.environ.get(name, "1").strip().lower() not in {"0", "false", "no", "off"}


def cache_path(home: Path | str | None = None) -> Path:
    """Location of the PubChem response cache inside the data directory."""
    return resolve_home(home) / "cache" / "pubchem-http.sqlite"


def cache_max_bytes() -> int:
    """The response cache's size limit in bytes (default 256 MB), set in megabytes
    with ``PUBCHEM_FINDER_CACHE_MAX_MB``. An unreadable value falls back to the
    default."""
    from .data.cache import DEFAULT_MAX_BYTES

    raw = os.environ.get(CACHE_MAX_MB_ENV)
    if raw is None:
        return DEFAULT_MAX_BYTES
    try:
        return max(0, round(float(raw) * 1_000_000))
    except ValueError:
        return DEFAULT_MAX_BYTES
