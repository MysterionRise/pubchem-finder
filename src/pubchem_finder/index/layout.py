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
import unicodedata
from collections.abc import Collection, Iterable
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, NamedTuple

from pydantic import BaseModel, Field

from .. import __version__
from ..core.exceptions import IndexIncompatible

FORMAT_VERSION = 1
"""Version of the index format this release writes and the newest it reads.

Bump it when an older pubchem-finder would misread an index written by this one
(a new optional field or Parquet column does not need it). Manifests without the
field were written by 2.1 or earlier: format 1.
"""


class SourceFile(BaseModel):
    """The local SDF file a chunk was built from, so a changed file is noticed."""

    path: str = Field(description="Resolved path")
    size: int
    mtime_ns: int

    @classmethod
    def of(cls, path: Path) -> "SourceFile":
        resolved = path.resolve()
        stat = resolved.stat()
        return cls(path=str(resolved), size=stat.st_size, mtime_ns=stat.st_mtime_ns)


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
    rdkit_version: str | None = Field(
        default=None,
        description="RDKit that fingerprinted the chunk (None: built by 2.1.0 or "
        "earlier; the chunk's .h5 records it)",
    )
    fpsim2_version: str | None = Field(
        default=None, description="FPSim2 that wrote the chunk's .h5"
    )
    source: SourceFile | None = Field(
        default=None,
        description="The file given with --sdf (None: a PubChem chunk, or a chunk "
        "built by 2.2.0 or earlier)",
    )
    skipped_id: str | None = Field(
        default=None,
        description="The first record ID that could not be used (IDs must be whole "
        "numbers from 1 to 4,294,967,295)",
    )
    build_id: str | None = Field(
        default=None,
        description="Unique to each build of the chunk; the manifest records which "
        "builds its merged files hold (None: built by 2.2.1 or earlier)",
    )


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
    fpsim2_version: str | None = Field(
        default=None, description="FPSim2 that wrote the .h5 (None: 2.1.0 or earlier)"
    )


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

    format_version: int = Field(
        default=1, description="Index format; see FORMAT_VERSION"
    )
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
    chunk_builds: dict[str, str | None] | None = Field(
        default=None,
        description="Each merged chunk's build_id, so a chunk rebuilt under its old "
        "name is noticed (None: merged by 2.2.1 or earlier)",
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
            ChunkStats.model_validate_json(p.read_text(encoding="utf-8"))
            for p in sorted(self.chunks.glob("*.json"))
        ]

    def completed_pattern_chunks(self) -> dict[str, PatternChunkStats]:
        if not self.pattern_chunks.is_dir():
            return {}
        stats = (
            PatternChunkStats.model_validate_json(p.read_text(encoding="utf-8"))
            for p in sorted(self.pattern_chunks.glob("*.json"))
        )
        return {s.name: s for s in stats}

    def load_manifest(self) -> Manifest | None:
        """The finalized index's manifest, or None if there is no index.

        Raises:
            IndexIncompatible: a newer pubchem-finder wrote the index.
        """
        merged = self.read_manifest()
        return merged.manifest if merged is not None else None

    def read_manifest(self) -> "Merged | None":
        """:meth:`load_manifest`, with the time the manifest was written."""
        if not self.fingerprints.exists():
            return None
        try:
            with self.manifest.open("rb") as stream:
                # Of the file being read, not of one a merge puts there meanwhile.
                written_ns = os.fstat(stream.fileno()).st_mtime_ns
                data = json.load(stream)
        except FileNotFoundError:  # no index, or a merge is swapping it right now
            return None
        # Before validating: a newer format may not fit this version's model.
        found = data.get("format_version", 1) if isinstance(data, dict) else 1
        if isinstance(found, int) and found > FORMAT_VERSION:
            built_by = data.get("pubchem_finder_version", "a newer version")
            raise IndexIncompatible(
                f"The index in {self.home} was built by pubchem-finder {built_by} "
                f"(index format {found}); pubchem-finder {__version__} reads format "
                f"{FORMAT_VERSION} and older. Upgrade pubchem-finder, or build an "
                "index for this version in another data directory (--home)."
            )
        return Merged(Manifest.model_validate(data), written_ns)

    def built_after(self, names: Iterable[str], when_ns: int) -> set[str]:
        """The chunks among ``names`` whose stats file was written after ``when_ns``."""
        out: set[str] = set()
        for name in names:
            try:
                if self.chunk_stats(name).stat().st_mtime_ns > when_ns:
                    out.add(name)
            except OSError:  # no stats: not a complete chunk
                continue
        return out


class Merged(NamedTuple):
    """A manifest as it was read, and when its file was written.

    A search keeps the one its fingerprints were loaded under and looks structures
    up against it: a merge that completes in between has put other records under
    the same IDs.
    """

    manifest: Manifest
    written_ns: int


def unmerged(
    manifest: Manifest, chunks: Iterable[ChunkStats], newer: Collection[str] = ()
) -> list[str]:
    """The chunks of a merged index that are no longer the builds it holds.

    ``chunks`` are the complete chunks on disk. A chunk rebuilt since the merge has
    another build id: its Parquet file describes the new file while the merged
    fingerprints are still the old one's, so the two must not be read together. A
    chunk without stats is being replaced right now, or is gone.

    ``newer`` names the chunks whose stats were written after the manifest (see
    :meth:`IndexPaths.built_after`). It decides for an index merged by 2.2.1 or
    earlier, which recorded no build ids.
    """
    on_disk = {c.name: c for c in chunks}
    out: list[str] = []
    for name in manifest.chunks:
        stats = on_disk.get(name)
        if stats is None:
            out.append(name)
        elif manifest.chunk_builds is None:
            # Merged by 2.2.1 or earlier. Those wrote the manifest last, so a chunk
            # of the user's own file with newer stats was built after the merge
            # (2.2.1 could rebuild one and fail before merging it). A build id
            # alone does not say so: 2.2.1 may have merged a chunk that has one.
            if stats.source is not None and name in newer:
                out.append(name)
        elif manifest.chunk_builds.get(name) != stats.build_id:
            out.append(name)
    return out


def chunk_key(name: str) -> str:
    """How chunk names are compared: ``Library`` and ``library`` are one set of files
    on macOS and Windows, so they are one chunk everywhere."""
    return unicodedata.normalize("NFC", name).casefold()


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
    tmp.write_text(text, encoding="utf-8")
    os.replace(tmp, path)
