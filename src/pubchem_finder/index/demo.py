"""The built-in demo index: FDA-approved drugs, built from data shipped in the package.

Local similarity and substructure search work on first use with it, before anyone
has built a PubChem index. It is a complete index of its own in
``<home>/demo/<key>/``, never mixed into the main index. The key hashes the data,
RDKit, FPSim2, the index format and the fingerprint settings, so an upgrade builds
a fresh demo index (in about a second) instead of reading a stale one.
"""

import contextlib
import csv
import hashlib
import io
import json
import logging
import os
import shutil
import sys
import threading
import uuid
from collections.abc import Iterator
from importlib.resources import files
from pathlib import Path
from typing import Any

from rdkit import Chem, rdBase

from ..config import data_dir_error
from ..core.exceptions import IndexBuildError, PubChemFinderError
from ..data.sdf import SdfRecord
from .builder import (
    DEFAULT_FP_PARAMS,
    FP_TYPE,
    EventFn,
    finalize,
    runtime_versions,
    write_chunk,
)
from .layout import FORMAT_VERSION, IndexPaths, Manifest

logger = logging.getLogger(__name__)

NAME = "fda_approved_drugs"
_DATA = files(__package__) / "demo_data"
_lock = threading.Lock()


def data() -> bytes:
    """The demo set: a TSV in the index's metadata columns."""
    return (_DATA / f"{NAME}.tsv").read_bytes()


def provenance() -> dict[str, Any]:
    """Where the demo set comes from (source, license, date, size)."""
    info: dict[str, Any] = json.loads((_DATA / f"{NAME}.json").read_text("utf-8"))
    return info


def key() -> str:
    """Changes whenever a demo index built now would differ from an older one."""
    digest = hashlib.sha256(data())
    settings = [*runtime_versions(), FORMAT_VERSION, FP_TYPE, DEFAULT_FP_PARAMS]
    digest.update(json.dumps(settings, sort_keys=True).encode())
    return digest.hexdigest()[:16]


def demo_home(home: Path) -> Path:
    """Where the demo index for this release lives inside the data directory."""
    return home / "demo" / key()


def manifest(home: Path) -> Manifest | None:
    """The demo index's manifest, or None if it isn't built yet."""
    return IndexPaths(demo_home(home)).load_manifest()


class TsvReader:
    """The demo TSV's rows as records for :func:`~.builder.write_chunk`."""

    def __init__(self, text: str) -> None:
        self.text = text
        self.n_ok = 0
        self.n_failed = 0

    def __iter__(self) -> Iterator[SdfRecord]:
        for row in csv.DictReader(io.StringIO(self.text), delimiter="\t"):
            mol = Chem.MolFromSmiles(row["smiles"]) if row["smiles"] else None
            if mol is None:
                self.n_failed += 1
                continue
            self.n_ok += 1
            yield SdfRecord(
                cid=int(row["cid"]),
                mol=mol,
                smiles=row["smiles"],
                inchikey=row["inchikey"] or None,
                iupac_name=row["iupac_name"] or None,
                molecular_formula=row["molecular_formula"] or None,
                molecular_weight=float(row["molecular_weight"])
                if row["molecular_weight"]
                else None,
            )


def ensure_built(home: Path, on_event: EventFn | None = None) -> Path:
    """The demo index's home, built first if needed.

    Builds in a temporary directory and renames it into place, so readers never
    see a half-built index; when another process wins the race, its index is used.
    Threads in one process build it once.
    """
    target = demo_home(home)
    if IndexPaths(target).load_manifest() is not None:
        return target
    with _lock:
        if IndexPaths(target).load_manifest() is not None:
            return target
        try:
            target.parent.mkdir(parents=True, exist_ok=True)
        except OSError as exc:
            raise data_dir_error(exc, home) from exc
        building = target.parent / f".building-{uuid.uuid4().hex}"
        try:
            _build(building, on_event)
            if target.exists():  # left by a build that was interrupted while moving in
                shutil.rmtree(target)
            try:
                os.rename(building, target)
            except OSError:
                if IndexPaths(target).load_manifest() is None:
                    raise
        except PubChemFinderError:
            raise
        except Exception as exc:
            raise IndexBuildError(f"Building the demo index failed: {exc}") from exc
        finally:
            shutil.rmtree(building, ignore_errors=True)
    return target


def _build(home: Path, on_event: EventFn | None) -> None:
    rdBase.DisableLog("rdApp.*")
    paths = IndexPaths(home)
    emit = on_event or (lambda msg: logger.info(msg))
    emit(f"building the demo index of {provenance()['name']}")
    # Over stdio, stdout is the MCP transport: nothing in the build may print to it.
    with contextlib.redirect_stdout(sys.stderr):
        reader = TsvReader(data().decode("utf-8"))
        write_chunk(paths, NAME, reader)
        finalize(paths, emit, substructure=True)
