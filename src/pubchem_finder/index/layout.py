"""On-disk layout of the local index and its manifest.

```
<home>/
  downloads/                 PubChem SDF chunks while they are being processed
  chunks/
    <chunk>.h5               FPSim2 fingerprints for one SDF chunk (unsorted)
    <chunk>.parquet          CID, SMILES, InChIKey, IUPAC name, formula, MW
    <chunk>.json             ChunkStats; written last, marks the chunk complete
    pattern/                 only with substructure search enabled
      <chunk>.h5             RDKit pattern fingerprints, built from <chunk>.parquet
      <chunk>.json           PatternChunkStats; written last
  index/
    pubchem.h5               merged, popcount-sorted FPSim2 database
    pubchem-pattern.h5       merged pattern fingerprints (substructure screening)
    manifest.json            Manifest; written last, marks the index complete
```
"""

import json
import os
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

from pydantic import BaseModel, Field


class ChunkStats(BaseModel):
    """Result of fingerprinting one SDF file."""

    name: str
    n_records: int = Field(description="Molecules fingerprinted and indexed")
    n_failed: int = Field(
        description="Records RDKit could not parse or that had no CID"
    )
    seconds: float
    fp_type: str
    fp_params: dict[str, Any]


def fpsim2_memory_bytes(n_compounds: int, fp_bits: int) -> int:
    """RAM FPSim2 needs to hold ``n_compounds`` fingerprints in memory
    (id + popcount + the fingerprint, as uint64 words)."""
    words = (fp_bits + 63) // 64 + 2
    return n_compounds * words * 8


PATTERN_FP_TYPE = "RDKitPattern"
PATTERN_FP_BITS = 2048  # FPSim2 always builds RDKit pattern fingerprints at 2048 bits


class PatternChunkStats(BaseModel):
    """Result of building pattern fingerprints for one chunk."""

    name: str
    n_records: int = Field(description="Molecules fingerprinted")
    n_failed: int = Field(description="Stored SMILES RDKit could not parse")
    seconds: float
    rdkit_version: str


class SubstructureInfo(BaseModel):
    """Describes the pattern-fingerprint database used for substructure search."""

    fp_type: str = PATTERN_FP_TYPE
    n_compounds: int
    n_failed: int
    chunks: list[str]
    built_at: datetime
    rdkit_version: str = Field(
        description="Pattern fingerprints can change between RDKit releases"
    )

    def memory_bytes(self) -> int:
        return fpsim2_memory_bytes(self.n_compounds, PATTERN_FP_BITS)


class Manifest(BaseModel):
    """Describes a finalized index."""

    fp_type: str
    fp_params: dict[str, Any]
    n_compounds: int
    n_failed: int
    chunks: list[str]
    built_at: datetime
    pubchem_finder_version: str
    rdkit_version: str
    fpsim2_version: str
    substructure: SubstructureInfo | None = Field(
        default=None, description="Present when substructure search is enabled"
    )

    @property
    def fingerprint(self) -> str:
        params = " ".join(f"{k}={v}" for k, v in sorted(self.fp_params.items()))
        return f"{self.fp_type} {params}".strip()

    def memory_bytes(self) -> int:
        """RAM needed to hold the similarity fingerprints in memory."""
        return fpsim2_memory_bytes(
            self.n_compounds, int(self.fp_params.get("fpSize", 2048))
        )


@dataclass(frozen=True)
class IndexPaths:
    """Paths inside a pubchem-finder data directory."""

    home: Path

    @property
    def downloads(self) -> Path:
        return self.home / "downloads"

    @property
    def chunks(self) -> Path:
        return self.home / "chunks"

    @property
    def index_dir(self) -> Path:
        return self.home / "index"

    @property
    def fingerprints(self) -> Path:
        return self.index_dir / "pubchem.h5"

    @property
    def pattern_fingerprints(self) -> Path:
        return self.index_dir / "pubchem-pattern.h5"

    @property
    def manifest(self) -> Path:
        return self.index_dir / "manifest.json"

    @property
    def pattern_chunks(self) -> Path:
        return self.chunks / "pattern"

    def chunk_h5(self, name: str) -> Path:
        return self.chunks / f"{name}.h5"

    def chunk_parquet(self, name: str) -> Path:
        return self.chunks / f"{name}.parquet"

    def chunk_stats(self, name: str) -> Path:
        return self.chunks / f"{name}.json"

    def pattern_h5(self, name: str) -> Path:
        return self.pattern_chunks / f"{name}.h5"

    def pattern_stats(self, name: str) -> Path:
        return self.pattern_chunks / f"{name}.json"

    def metadata_files(self, names: list[str]) -> list[Path]:
        return [p for p in (self.chunk_parquet(n) for n in names) if p.exists()]

    def completed_chunks(self) -> list[ChunkStats]:
        if not self.chunks.is_dir():
            return []
        return [
            ChunkStats.model_validate_json(p.read_text())
            for p in sorted(self.chunks.glob("*.json"))
        ]

    def completed_pattern_chunks(self) -> dict[str, PatternChunkStats]:
        if not self.pattern_chunks.is_dir():
            return {}
        stats = (
            PatternChunkStats.model_validate_json(p.read_text())
            for p in sorted(self.pattern_chunks.glob("*.json"))
        )
        return {s.name: s for s in stats}

    def load_manifest(self) -> Manifest | None:
        if not self.manifest.exists() or not self.fingerprints.exists():
            return None
        return Manifest.model_validate_json(self.manifest.read_text())


def chunk_name(path: Path | str) -> str:
    """``Compound_000000001_000500000.sdf.gz`` -> ``Compound_000000001_000500000``."""
    name = Path(path).name
    for suffix in (".gz", ".sdf", ".sd", ".mol"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
    return name


def write_json_atomic(path: Path, data: BaseModel | dict[str, Any]) -> None:
    """Write JSON to ``path`` via a temporary file and rename."""
    text = (
        data.model_dump_json(indent=2)
        if isinstance(data, BaseModel)
        else json.dumps(data, indent=2)
    )
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(text)
    os.replace(tmp, path)
