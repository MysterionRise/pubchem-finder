"""Pydantic models shared by the CLI and the MCP server."""

from pydantic import BaseModel, Field, computed_field

PUBCHEM_COMPOUND_URL = "https://pubchem.ncbi.nlm.nih.gov/compound/{cid}"


class CompoundRecord(BaseModel):
    """One PubChem compound as returned by the PUG REST property table."""

    cid: int = Field(description="PubChem Compound ID")
    title: str | None = Field(None, description="PubChem title (common name)")
    iupac_name: str | None = None
    smiles: str | None = Field(None, description="SMILES with stereo and isotopes")
    connectivity_smiles: str | None = Field(
        None, description="SMILES without stereo/isotopes"
    )
    inchikey: str | None = None
    molecular_formula: str | None = None
    molecular_weight: float | None = None
    xlogp: float | None = None
    h_bond_donors: int | None = None
    h_bond_acceptors: int | None = None
    tpsa: float | None = None

    @computed_field  # type: ignore[prop-decorator]
    @property
    def pubchem_url(self) -> str:
        return PUBCHEM_COMPOUND_URL.format(cid=self.cid)


class LookupResult(BaseModel):
    """Result of a PubChem lookup."""

    query: str
    by: str
    compounds: list[CompoundRecord]


class Properties(BaseModel):
    """RDKit descriptors for one structure."""

    input_smiles: str
    canonical_smiles: str
    molecular_weight: float
    logp: float
    num_h_donors: int
    num_h_acceptors: int
    tpsa: float
    num_rotatable_bonds: int
    num_aromatic_rings: int
    num_heavy_atoms: int
    lipinski_ok: bool = Field(description="Passes Lipinski's Rule of Five")


class Standardized(BaseModel):
    """Result of standardizing a structure."""

    input_smiles: str
    smiles: str = Field(description="Canonical SMILES of the standardized parent")
    changed: bool = Field(description="Whether standardization altered the structure")


class Comparison(BaseModel):
    """Fingerprint similarity between two structures."""

    smiles_a: str
    smiles_b: str
    ecfp4_tanimoto: float = Field(description="Tanimoto on Morgan r=2, 2048 bits")
    maccs_tanimoto: float = Field(description="Tanimoto on MACCS keys")


class SimilarityHit(BaseModel):
    """One compound from the local index."""

    cid: int
    similarity: float = Field(description="Similarity to the query, 0 to 1")
    smiles: str | None = None
    inchikey: str | None = None
    iupac_name: str | None = None
    molecular_formula: str | None = None
    molecular_weight: float | None = None

    @computed_field  # type: ignore[prop-decorator]
    @property
    def pubchem_url(self) -> str:
        return PUBCHEM_COMPOUND_URL.format(cid=self.cid)


class SimilarityResults(BaseModel):
    """Result of a similarity search over the local index."""

    query_smiles: str
    threshold: float
    metric: str
    fingerprint: str = Field(description="Fingerprint the index was built with")
    n_indexed: int = Field(description="Compounds in the index")
    hits: list[SimilarityHit]


class IndexStatus(BaseModel):
    """State of the local index."""

    state: str = Field(description="'ready' or 'not_built'")
    home: str = Field(description="Data directory")
    n_compounds: int = 0
    n_failed: int = 0
    fingerprint: str | None = None
    built_at: str | None = None
    size_bytes: int = 0
    memory_bytes: int = Field(default=0, description="RAM needed to search in memory")
    chunks_in_index: int = 0
    chunks_built: int = Field(default=0, description="Chunks on disk, merged or not")
    loaded: bool = Field(
        default=False, description="Whether this process has the index in use"
    )
    mode: str | None = Field(default=None, description="'memory' or 'disk' once loaded")
    hint: str | None = None
