"""PubChem Finder: molecular similarity search over PubChem."""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("pubchem-finder")
except PackageNotFoundError:  # pragma: no cover - source checkout without install
    __version__ = "0.0.0"

__all__ = ["__version__"]
