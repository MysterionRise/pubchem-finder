"""On-disk layout of the local index and its manifest.

```
<home>/
  downloads/                 PubChem SDF chunks while they are being processed
  chunks/
    <chunk>.h5               FPSim2 fingerprints for one SDF chunk (unsorted)
    <chunk>.parquet          CID, SMILES, InChIKey, IUPAC name, formula, MW
    <chunk>.json             ChunkStats; written last, marks the chunk complete
  index/
    pubchem.h5               merged, popcount-sorted FPSim2 database
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

    @property
    def fingerprint(self) -> str:
        params = " ".join(f"{k}={v}" for k, v in sorted(self.fp_params.items()))
        return f"{self.fp_type} {params}".strip()

    def memory_bytes(self) -> int:
        """RAM needed to hold the fingerprints in memory (FPSim2 layout:
        id + popcount + the fingerprint, as uint64 words)."""
        fp_bits = int(self.fp_params.get("fpSize", 2048))
        words = (fp_bits + 63) // 64 + 2
        return self.n_compounds * words * 8


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
    def manifest(self) -> Path:
        return self.index_dir / "manifest.json"

    def chunk_h5(self, name: str) -> Path:
        return self.chunks / f"{name}.h5"

    def chunk_parquet(self, name: str) -> Path:
        return self.chunks / f"{name}.parquet"

    def chunk_stats(self, name: str) -> Path:
        return self.chunks / f"{name}.json"

    def metadata_files(self, names: list[str]) -> list[Path]:
        return [p for p in (self.chunk_parquet(n) for n in names) if p.exists()]

    def completed_chunks(self) -> list[ChunkStats]:
        if not self.chunks.is_dir():
            return []
        return [
            ChunkStats.model_validate_json(p.read_text())
            for p in sorted(self.chunks.glob("*.json"))
        ]

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
