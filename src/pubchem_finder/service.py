"""Framework-free operations shared by the CLI and the MCP server.

Every function here takes plain values, returns a Pydantic model and raises
``PubChemFinderError`` subclasses. The CLI and the MCP server only translate
those into their own output and error formats.
"""

import threading
from pathlib import Path

from rdkit import Chem, DataStructs, rdBase

from .chemistry import MoleculeStandardizer, calculate_properties, is_druglike
from .config import resolve_home
from .core.exceptions import MoleculeValidationError, PubChemFinderError
from .data.pug_rest import Namespace, PubChemClient
from .fingerprints import MACCSFingerprint, MorganFingerprint
from .index.engine import Metric, SimilarityIndex
from .models import (
    Comparison,
    IndexStatus,
    LookupResult,
    Properties,
    SimilarityHit,
    SimilarityResults,
    Standardized,
)

# RDKit writes parse warnings to stderr; keep server output quiet.
rdBase.DisableLog("rdApp.*")

_ECFP4 = MorganFingerprint(radius=2, n_bits=2048)
_MACCS = MACCSFingerprint()
_STANDARDIZER = MoleculeStandardizer()

_client_lock = threading.Lock()
_client: PubChemClient | None = None
_indexes: dict[Path, SimilarityIndex] = {}


def get_client() -> PubChemClient:
    """Process-wide PubChem client, created on first use."""
    global _client
    with _client_lock:
        if _client is None:
            _client = PubChemClient()
        return _client


def parse_smiles(smiles: str) -> Chem.Mol:
    """Parse SMILES or raise ``MoleculeValidationError``."""
    mol: Chem.Mol | None = (
        Chem.MolFromSmiles(smiles.strip()) if smiles.strip() else None
    )
    if mol is None:
        raise MoleculeValidationError(f"Invalid SMILES: {smiles!r}")
    return mol


def lookup(query: str, by: Namespace = "name", limit: int = 5) -> LookupResult:
    """Look compounds up in PubChem."""
    compounds = get_client().lookup(query, by=by, limit=limit)
    return LookupResult(query=query, by=by, compounds=compounds)


def properties(smiles: str) -> Properties:
    """RDKit descriptors and Lipinski check for one structure."""
    mol = parse_smiles(smiles)
    props = calculate_properties(mol)
    return Properties.model_validate(
        {
            "input_smiles": smiles,
            "canonical_smiles": Chem.MolToSmiles(mol),
            "lipinski_ok": is_druglike(mol),
            **props,
        }
    )


def standardize(smiles: str) -> Standardized:
    """Standardize a structure (cleanup, parent fragment, neutralize)."""
    mol = parse_smiles(smiles)
    std = _STANDARDIZER.standardize(mol)
    if std is None:
        raise MoleculeValidationError(f"Could not standardize {smiles!r}")
    result = str(Chem.MolToSmiles(std))
    return Standardized(
        input_smiles=smiles,
        smiles=result,
        changed=result != Chem.MolToSmiles(mol),
    )


def compare(smiles_a: str, smiles_b: str) -> Comparison:
    """Tanimoto similarity between two structures on ECFP4 and MACCS."""
    mol_a = parse_smiles(smiles_a)
    mol_b = parse_smiles(smiles_b)
    return Comparison(
        smiles_a=smiles_a,
        smiles_b=smiles_b,
        ecfp4_tanimoto=float(
            DataStructs.TanimotoSimilarity(
                _ECFP4.generate_bitvect(mol_a), _ECFP4.generate_bitvect(mol_b)
            )
        ),
        maccs_tanimoto=float(
            DataStructs.TanimotoSimilarity(
                _MACCS.generate_bitvect(mol_a), _MACCS.generate_bitvect(mol_b)
            )
        ),
    )


def get_index(home: Path | str | None = None) -> SimilarityIndex:
    """The local index for ``home`` (default: configured data directory), cached."""
    path = resolve_home(home)
    with _client_lock:
        if path not in _indexes:
            _indexes[path] = SimilarityIndex(path)
        return _indexes[path]


def similarity_search(
    smiles: str,
    threshold: float = 0.7,
    top_k: int = 20,
    metric: Metric = "tanimoto",
    home: Path | str | None = None,
) -> SimilarityResults:
    """Most similar compounds in the local index."""
    if not 0.0 < threshold <= 1.0:
        raise PubChemFinderError("threshold must be in (0, 1]")
    if not 1 <= top_k <= 1000:
        raise PubChemFinderError("top_k must be between 1 and 1000")
    if metric not in ("tanimoto", "dice", "cosine"):
        raise PubChemFinderError(f"Unknown metric {metric!r}")
    mol = parse_smiles(smiles)
    index = get_index(home)
    manifest = index.manifest
    pairs = index.search(mol, threshold=threshold, top_k=top_k, metric=metric)
    meta = index.metadata([cid for cid, _ in pairs])
    hits = [
        SimilarityHit.model_validate(
            {**meta.get(cid, {}), "cid": cid, "similarity": sim}
        )
        for cid, sim in pairs
    ]
    return SimilarityResults(
        query_smiles=smiles,
        threshold=threshold,
        metric=metric,
        fingerprint=manifest.fingerprint,
        n_indexed=manifest.n_compounds,
        hits=hits,
    )


def index_status(home: Path | str | None = None) -> IndexStatus:
    """Describe the local index without loading it."""
    index = get_index(home)
    paths = index.paths
    built = paths.completed_chunks()
    manifest = paths.load_manifest()
    if manifest is None:
        return IndexStatus(
            state="not_built",
            home=str(paths.home),
            chunks_built=len(built),
            hint=(
                "Build with `pubchem-finder index build --limit 1` for a ~500k-compound "
                "sample, or `--all` for all of PubChem."
            ),
        )
    return IndexStatus(
        state="ready",
        home=str(paths.home),
        n_compounds=manifest.n_compounds,
        n_failed=manifest.n_failed,
        fingerprint=manifest.fingerprint,
        built_at=manifest.built_at.isoformat(timespec="seconds"),
        size_bytes=paths.fingerprints.stat().st_size,
        memory_bytes=manifest.memory_bytes(),
        chunks_in_index=len(manifest.chunks),
        chunks_built=len(built),
        loaded=index.loaded,
        mode=None
        if index.in_memory is None
        else ("memory" if index.in_memory else "disk"),
    )
