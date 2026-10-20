"""Local similarity index: build from SDF files, search with FPSim2."""

from .builder import BuildReport, SdfSource, build_chunks, finalize
from .engine import SimilarityIndex
from .layout import ChunkStats, IndexPaths, Manifest

__all__ = [
    "BuildReport",
    "ChunkStats",
    "IndexPaths",
    "Manifest",
    "SdfSource",
    "SimilarityIndex",
    "build_chunks",
    "finalize",
]
