"""Framework-free operations shared by the CLI and the MCP server.

Every function here takes plain values, returns a Pydantic model and raises
``PubChemFinderError`` subclasses. The CLI and the MCP server only translate
those into their own output and error formats.
"""

import threading

from rdkit import Chem, DataStructs, rdBase

from .chemistry import MoleculeStandardizer, calculate_properties, is_druglike
from .core.exceptions import MoleculeValidationError
from .data.pug_rest import Namespace, PubChemClient
from .fingerprints import MACCSFingerprint, MorganFingerprint
from .models import Comparison, LookupResult, Properties, Standardized

# RDKit writes parse warnings to stderr; keep server output quiet.
rdBase.DisableLog("rdApp.*")

_ECFP4 = MorganFingerprint(radius=2, n_bits=2048)
_MACCS = MACCSFingerprint()
_STANDARDIZER = MoleculeStandardizer()

_client_lock = threading.Lock()
_client: PubChemClient | None = None


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
