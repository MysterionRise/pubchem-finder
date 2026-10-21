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
    scaffolds: list["ScaffoldGroup"] | None = Field(
        default=None, description="Hits grouped by Murcko scaffold, when requested"
    )


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


# --------------------------------------------------------------------------- #
# PubChem knowledge tools (Phase 4a)
# --------------------------------------------------------------------------- #


class RemoteHit(CompoundRecord):
    """A compound found by one of PubChem's structure searches."""

    ecfp4_similarity: float | None = Field(
        default=None,
        description=(
            "Tanimoto similarity to the query on ECFP4 (Morgan radius 2, 2048 bits), "
            "computed locally; comparable with the local similarity_search"
        ),
    )


class RemoteSearchResults(BaseModel):
    """Result of a PubChem structure search over all of PubChem."""

    query_smiles: str
    search: str = Field(description="'similarity', 'substructure' or 'identity'")
    threshold: int | None = Field(
        default=None, description="PubChem similarity threshold, 0-100"
    )
    max_records: int
    truncated: bool = Field(
        description="PubChem stopped at max_records; more may exist"
    )
    hits: list[RemoteHit]
    note: str


class Synonyms(BaseModel):
    """Names PubChem knows a compound by, most common first."""

    cid: int
    synonyms: list[str]
    cas_numbers: list[str] = Field(
        description="CAS Registry Numbers with a valid check digit"
    )
    total: int = Field(description="Synonyms PubChem has before the limit was applied")


class CrossRefs(BaseModel):
    """External identifiers linked to a compound."""

    cid: int
    counts: dict[str, int] = Field(description="Number of identifiers of each type")
    samples: dict[str, list[str]] = Field(
        description="First few identifiers of each type"
    )
    pubmed_urls: list[str] = Field(default_factory=list)


class SourcedText(BaseModel):
    text: str
    source: str | None = None
    url: str | None = None


class Description(BaseModel):
    """Short curated descriptions of a compound."""

    cid: int
    title: str | None
    descriptions: list[SourcedText]


class TargetActivity(BaseModel):
    target_accession: str | None
    target_gene_id: str | None
    active_assays: int


class PotentResult(BaseModel):
    aid: int
    assay_name: str | None
    activity_name: str | None
    value_um: float
    target_accession: str | None


class BioactivitySummary(BaseModel):
    """Aggregate of PubChem's bioassay results for one compound."""

    cid: int
    tested: bool
    n_assays: int = Field(description="Distinct BioAssay AIDs")
    outcomes: dict[str, int] = Field(description="Result counts by outcome")
    top_targets: list[TargetActivity] = Field(
        description="Targets with the most active results"
    )
    most_potent: list[PotentResult] = Field(
        description="Active results with the lowest activity values (µM)"
    )

    @computed_field  # type: ignore[prop-decorator]
    @property
    def pubchem_url(self) -> str:
        return PUBCHEM_COMPOUND_URL.format(cid=self.cid) + "#section=BioAssay-Results"


class Annotation(BaseModel):
    section: str = Field(
        description="Section path, e.g. 'Hazards Identification > GHS Classification'"
    )
    name: str | None = None
    text: str
    source: str | None = None
    url: str | None = None
    license: str | None = None


class AnnotationSection(BaseModel):
    """Curated annotations from a PubChem compound page section."""

    cid: int
    title: str | None
    heading: str
    annotations: list[Annotation]
    n_total: int = Field(description="Annotations available before trimming")
    truncated: bool
    sources: list[str] = Field(description="Depositors of the annotations shown")


# --------------------------------------------------------------------------- #
# Assistant features (Phase 4b)
# --------------------------------------------------------------------------- #


class StructuralAlert(BaseModel):
    """A substructure flagged by a published filter set."""

    catalog: str = Field(description="Filter set, e.g. PAINS_A, Brenk, NIH")
    name: str = Field(description="Alert name, e.g. ene_rhod_A(235)")
    atoms: list[int] = Field(
        description="Indices of the matched atoms in the input SMILES"
    )
    reference: str | None = Field(
        default=None, description="Publication defining the filter"
    )


class MedchemProfile(BaseModel):
    """Drug-likeness and liability triage for one structure."""

    input_smiles: str
    canonical_smiles: str
    qed: float = Field(
        description="Quantitative estimate of drug-likeness, 0 (poor) to 1"
    )
    sa_score: float | None = Field(
        default=None,
        description="Synthetic accessibility, 1 (easy) to 10 (very hard); None if unavailable",
    )
    molecular_weight: float
    logp: float
    num_h_donors: int
    num_h_acceptors: int
    tpsa: float
    num_rotatable_bonds: int
    lipinski_violations: list[str] = Field(description="Broken Rule-of-Five criteria")
    veber_ok: bool = Field(description="Rotatable bonds <= 10 and TPSA <= 140")
    alerts: list[StructuralAlert]
    pains_free: bool = Field(description="No PAINS alert matched")
    murcko_scaffold: str = Field(
        description="Bemis-Murcko scaffold SMILES; empty if acyclic"
    )


class ScaffoldGroup(BaseModel):
    """Hits that share a Bemis-Murcko scaffold."""

    scaffold: str = Field(description="Scaffold SMILES; empty for acyclic compounds")
    cids: list[int]
    best_similarity: float
