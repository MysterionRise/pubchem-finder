"""Local similarity index: build from SDF files, search with FPSim2."""

from .builder import BuildReport, SdfSource, build_chunks, finalize
from .engine import LocalIndex
from .layout import ChunkStats, IndexPaths, Manifest

__all__ = [
    "BuildReport",
    "ChunkStats",
    "IndexPaths",
    "LocalIndex",
    "Manifest",
    "SdfSource",
    "build_chunks",
    "finalize",
]
