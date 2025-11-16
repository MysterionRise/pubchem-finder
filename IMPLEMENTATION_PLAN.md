# PubChem Finder Implementation Plan
## Building a Solid Foundation for Advanced Molecular Similarity Search

---

## Philosophy: Foundation First, Extension Later

This implementation plan follows the principle of **"Build the foundation right, extend incrementally"**:

1. **Solid Core**: Start with clean, well-tested core components
2. **Plugin Architecture**: Design for extensibility without modification
3. **Gradual Enhancement**: Add complexity only when needed
4. **Production Ready**: Each phase delivers working, deployable software

---

## Architecture Overview: Layered & Extensible

```
┌─────────────────────────────────────────────────────────────┐
│                     Future: ML Layer                        │
│          (GNNs, Transformers - Pluggable)                   │
└─────────────────────────────────────────────────────────────┘
                            ▲
                            │ Extension Point
┌─────────────────────────────────────────────────────────────┐
│              Phase 3: API & Advanced Search                 │
│         (FastAPI, Multiple Search Modes)                    │
└─────────────────────────────────────────────────────────────┘
                            ▲
                            │ Built Upon
┌─────────────────────────────────────────────────────────────┐
│          Phase 2: Basic Search & Vector Storage             │
│    (Simple Similarity, Vector DB Integration)               │
└─────────────────────────────────────────────────────────────┘
                            ▲
                            │ Built Upon
┌─────────────────────────────────────────────────────────────┐
│            Phase 1: CORE FOUNDATION (START HERE)            │
│  (RDKit, Fingerprints, Data Models, Ingestion)              │
└─────────────────────────────────────────────────────────────┘
```

---

## Phase 1: Core Foundation (Weeks 1-4) 🎯 **START HERE**

### Goal
Build robust, well-tested core components that will support all future features.

### Deliverables
✅ Clean project structure with proper packaging
✅ RDKit-based molecule processing pipeline
✅ Multiple fingerprint generators (extensible design)
✅ Relational database schema and models
✅ PubChem data ingestion pipeline
✅ Comprehensive test suite
✅ CLI for basic operations

---

### 1.1 Project Setup & Structure

#### Directory Structure
```
pubchem-finder/
├── pyproject.toml              # Modern Python packaging (Poetry/PDM)
├── README.md
├── CHANGELOG.md
├── docker-compose.yml          # Development environment
├── .env.example
│
├── src/
│   └── pubchem_finder/
│       ├── __init__.py
│       │
│       ├── core/                      # Core abstractions
│       │   ├── __init__.py
│       │   ├── molecule.py            # Molecule data class
│       │   ├── fingerprint.py         # Fingerprint base classes
│       │   └── exceptions.py          # Custom exceptions
│       │
│       ├── chemistry/                 # Chemistry processing
│       │   ├── __init__.py
│       │   ├── standardizer.py        # Molecule standardization
│       │   ├── validator.py           # Validation logic
│       │   └── properties.py          # Property calculation
│       │
│       ├── fingerprints/              # Fingerprint implementations
│       │   ├── __init__.py
│       │   ├── base.py                # Abstract base
│       │   ├── morgan.py              # Morgan/ECFP
│       │   ├── maccs.py               # MACCS keys
│       │   ├── rdkit_fp.py            # RDKit topological
│       │   └── registry.py            # Fingerprint registry
│       │
│       ├── database/                  # Data persistence
│       │   ├── __init__.py
│       │   ├── models.py              # SQLAlchemy models
│       │   ├── repository.py          # Data access layer
│       │   └── migrations/            # Alembic migrations
│       │
│       ├── ingestion/                 # Data ingestion
│       │   ├── __init__.py
│       │   ├── pubchem.py             # PubChem FTP downloader
│       │   ├── sdf_parser.py          # SDF file parser
│       │   └── pipeline.py            # Ingestion pipeline
│       │
│       ├── config/                    # Configuration
│       │   ├── __init__.py
│       │   └── settings.py            # Pydantic settings
│       │
│       └── cli/                       # Command-line interface
│           ├── __init__.py
│           └── commands.py            # Click commands
│
├── tests/
│   ├── unit/                          # Unit tests
│   ├── integration/                   # Integration tests
│   ├── fixtures/                      # Test data
│   └── conftest.py                    # Pytest configuration
│
├── scripts/                           # Utility scripts
│   ├── init_db.py
│   └── benchmark.py
│
└── docs/                              # Documentation
    ├── architecture.md
    └── api/
```

#### Task 1.1: Initial Setup
```bash
# What to implement:
1. Create new branch: feature/phase1-foundation
2. Set up pyproject.toml with Poetry or PDM
3. Configure development tools:
   - black (code formatting)
   - ruff (fast linting)
   - mypy (type checking)
   - pytest (testing)
4. Create directory structure
5. Set up pre-commit hooks
```

**Implementation Time**: 1 day

---

### 1.2 Core Data Models

#### Define Core Abstractions

**File: `src/pubchem_finder/core/molecule.py`**
```python
"""Core molecule data structures."""
from dataclasses import dataclass, field
from typing import Optional, Dict, Any
from datetime import datetime


@dataclass
class Molecule:
    """
    Core molecule data structure.

    This is the central data model used throughout the application.
    Intentionally simple and focused on essentials.
    """
    # Identifiers
    id: Optional[str] = None  # UUID or PubChem CID

    # Chemical representations
    smiles: str = ""
    canonical_smiles: str = ""
    inchi: Optional[str] = None
    inchi_key: Optional[str] = None

    # Molecular properties (calculated)
    molecular_weight: Optional[float] = None
    logp: Optional[float] = None
    num_h_donors: Optional[int] = None
    num_h_acceptors: Optional[int] = None
    tpsa: Optional[float] = None  # Topological polar surface area
    num_rotatable_bonds: Optional[int] = None
    num_aromatic_rings: Optional[int] = None

    # Fingerprints (stored as dict: {type: vector})
    fingerprints: Dict[str, Any] = field(default_factory=dict)

    # Metadata
    source: str = "pubchem"
    metadata: Dict[str, Any] = field(default_factory=dict)
    created_at: Optional[datetime] = None
    updated_at: Optional[datetime] = None

    def __post_init__(self):
        """Validate essential fields."""
        if not self.smiles and not self.canonical_smiles:
            raise ValueError("Either smiles or canonical_smiles must be provided")
```

**File: `src/pubchem_finder/core/fingerprint.py`**
```python
"""Abstract base classes for fingerprints."""
from abc import ABC, abstractmethod
from typing import Any
import numpy as np
from rdkit import Chem


class FingerprintGenerator(ABC):
    """
    Abstract base class for all fingerprint generators.

    This interface ensures all fingerprint implementations are interchangeable.
    Future fingerprint types (including ML-based) will implement this interface.
    """

    @property
    @abstractmethod
    def name(self) -> str:
        """Unique identifier for this fingerprint type."""
        pass

    @property
    @abstractmethod
    def size(self) -> int:
        """Dimensionality of the fingerprint."""
        pass

    @abstractmethod
    def generate(self, mol: Chem.Mol) -> np.ndarray:
        """
        Generate fingerprint from RDKit molecule.

        Args:
            mol: RDKit Mol object

        Returns:
            numpy array representing the fingerprint

        Raises:
            FingerprintGenerationError: If fingerprint cannot be generated
        """
        pass

    def generate_from_smiles(self, smiles: str) -> np.ndarray:
        """Convenience method to generate from SMILES."""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {smiles}")
        return self.generate(mol)


class SimilarityMetric(ABC):
    """
    Abstract base for similarity metrics.

    Extensible design allows adding new metrics without changing search logic.
    """

    @property
    @abstractmethod
    def name(self) -> str:
        """Metric identifier."""
        pass

    @abstractmethod
    def calculate(self, fp1: np.ndarray, fp2: np.ndarray) -> float:
        """
        Calculate similarity between two fingerprints.

        Returns:
            Similarity score (higher = more similar)
        """
        pass
```

**Implementation Time**: 2 days

---

### 1.3 Chemistry Processing Layer

#### Molecule Standardization

**File: `src/pubchem_finder/chemistry/standardizer.py`**
```python
"""Molecule standardization and normalization."""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize
from typing import Optional
import logging

logger = logging.getLogger(__name__)


class MoleculeStandardizer:
    """
    Standardize molecules to canonical form.

    This is CRITICAL for consistent fingerprints and searches.
    All molecules should be standardized before fingerprinting.
    """

    def __init__(self):
        self.normalizer = rdMolStandardize.Normalizer()
        self.uncharger = rdMolStandardize.Uncharger()

    def standardize(self, mol: Chem.Mol) -> Optional[Chem.Mol]:
        """
        Standardize a molecule.

        Steps:
        1. Remove hydrogen atoms
        2. Normalize functional groups
        3. Neutralize charges (optional)
        4. Select largest fragment
        5. Canonicalize
        """
        try:
            # Remove explicit hydrogens
            mol = Chem.RemoveHs(mol)

            # Normalize (standardize functional groups)
            mol = self.normalizer.normalize(mol)

            # Remove charges (makes similarity more robust)
            mol = self.uncharger.uncharge(mol)

            # Take largest fragment (removes salts, counterions)
            mol = self._get_largest_fragment(mol)

            # Canonicalize
            Chem.SanitizeMol(mol)

            return mol

        except Exception as e:
            logger.error(f"Standardization failed: {e}")
            return None

    def _get_largest_fragment(self, mol: Chem.Mol) -> Chem.Mol:
        """Select largest fragment from molecule."""
        frags = Chem.GetMolFrags(mol, asMols=True)
        if len(frags) == 1:
            return mol
        # Return fragment with most heavy atoms
        return max(frags, key=lambda m: m.GetNumHeavyAtoms())

    def standardize_smiles(self, smiles: str) -> Optional[str]:
        """Standardize SMILES string."""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        mol = self.standardize(mol)
        if mol is None:
            return None

        return Chem.MolToSmiles(mol)
```

#### Property Calculator

**File: `src/pubchem_finder/chemistry/properties.py`**
```python
"""Calculate molecular descriptors and properties."""
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, Crippen
from typing import Dict


def calculate_properties(mol: Chem.Mol) -> Dict[str, float]:
    """
    Calculate common molecular properties.

    These properties can be used for filtering search results.
    """
    return {
        "molecular_weight": Descriptors.MolWt(mol),
        "logp": Crippen.MolLogP(mol),
        "num_h_donors": Lipinski.NumHDonors(mol),
        "num_h_acceptors": Lipinski.NumHAcceptors(mol),
        "tpsa": Descriptors.TPSA(mol),
        "num_rotatable_bonds": Lipinski.NumRotatableBonds(mol),
        "num_aromatic_rings": Lipinski.NumAromaticRings(mol),
        "num_heavy_atoms": Lipinski.HeavyAtomCount(mol),
    }
```

**Implementation Time**: 2 days

---

### 1.4 Fingerprint Generators

#### Base Implementation

**File: `src/pubchem_finder/fingerprints/base.py`**
```python
"""Base fingerprint implementations."""
from ..core.fingerprint import FingerprintGenerator, SimilarityMetric
import numpy as np


class TanimotoSimilarity(SimilarityMetric):
    """Tanimoto coefficient for binary fingerprints."""

    @property
    def name(self) -> str:
        return "tanimoto"

    def calculate(self, fp1: np.ndarray, fp2: np.ndarray) -> float:
        """Calculate Tanimoto similarity."""
        intersection = np.sum(fp1 & fp2)
        union = np.sum(fp1 | fp2)
        return float(intersection / union) if union > 0 else 0.0


class DiceSimilarity(SimilarityMetric):
    """Dice coefficient."""

    @property
    def name(self) -> str:
        return "dice"

    def calculate(self, fp1: np.ndarray, fp2: np.ndarray) -> float:
        """Calculate Dice similarity."""
        intersection = np.sum(fp1 & fp2)
        total = np.sum(fp1) + np.sum(fp2)
        return float(2 * intersection / total) if total > 0 else 0.0
```

**File: `src/pubchem_finder/fingerprints/morgan.py`**
```python
"""Morgan (ECFP) fingerprint implementation."""
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
from ..core.fingerprint import FingerprintGenerator


class MorganFingerprint(FingerprintGenerator):
    """
    Morgan fingerprint generator (aka ECFP).

    radius=2 corresponds to ECFP4 (diameter=4)
    radius=3 corresponds to ECFP6 (diameter=6)
    """

    def __init__(self, radius: int = 2, n_bits: int = 2048, use_features: bool = False):
        self.radius = radius
        self.n_bits = n_bits
        self.use_features = use_features

    @property
    def name(self) -> str:
        feature_str = "fcfp" if self.use_features else "ecfp"
        return f"morgan_{feature_str}_r{self.radius}_{self.n_bits}"

    @property
    def size(self) -> int:
        return self.n_bits

    def generate(self, mol: Chem.Mol) -> np.ndarray:
        """Generate Morgan fingerprint."""
        fp = AllChem.GetMorganFingerprintAsBitVect(
            mol,
            radius=self.radius,
            nBits=self.n_bits,
            useFeatures=self.use_features
        )
        # Convert to numpy array
        arr = np.zeros((self.n_bits,), dtype=np.uint8)
        AllChem.DataStructs.ConvertToNumpyArray(fp, arr)
        return arr
```

**File: `src/pubchem_finder/fingerprints/maccs.py`**
```python
"""MACCS keys fingerprint implementation."""
from rdkit import Chem
from rdkit.Chem import MACCSkeys
import numpy as np
from ..core.fingerprint import FingerprintGenerator


class MACCSFingerprint(FingerprintGenerator):
    """MACCS 166-bit structural keys."""

    @property
    def name(self) -> str:
        return "maccs"

    @property
    def size(self) -> int:
        return 167  # MACCS keys are 167 bits (0-166, but 0 is not used)

    def generate(self, mol: Chem.Mol) -> np.ndarray:
        """Generate MACCS fingerprint."""
        fp = MACCSkeys.GenMACCSKeys(mol)
        arr = np.zeros((167,), dtype=np.uint8)
        for i in range(167):
            arr[i] = fp[i]
        return arr
```

**File: `src/pubchem_finder/fingerprints/registry.py`**
```python
"""Fingerprint registry for dynamic fingerprint selection."""
from typing import Dict, Type
from ..core.fingerprint import FingerprintGenerator
from .morgan import MorganFingerprint
from .maccs import MACCSFingerprint


class FingerprintRegistry:
    """
    Registry pattern for fingerprint generators.

    This makes it easy to add new fingerprint types without modifying
    existing code. Future ML-based fingerprints can be registered here.
    """

    def __init__(self):
        self._generators: Dict[str, FingerprintGenerator] = {}
        self._register_defaults()

    def _register_defaults(self):
        """Register standard fingerprints."""
        # Morgan variants
        self.register("morgan_2", MorganFingerprint(radius=2, n_bits=2048))
        self.register("morgan_3", MorganFingerprint(radius=3, n_bits=2048))
        self.register("ecfp4", MorganFingerprint(radius=2, n_bits=2048))  # Alias
        self.register("ecfp6", MorganFingerprint(radius=3, n_bits=2048))  # Alias

        # MACCS
        self.register("maccs", MACCSFingerprint())

    def register(self, name: str, generator: FingerprintGenerator):
        """Register a new fingerprint generator."""
        self._generators[name] = generator

    def get(self, name: str) -> FingerprintGenerator:
        """Get fingerprint generator by name."""
        if name not in self._generators:
            raise KeyError(f"Unknown fingerprint type: {name}")
        return self._generators[name]

    def list_available(self) -> list[str]:
        """List all registered fingerprint types."""
        return list(self._generators.keys())


# Global registry instance
fingerprint_registry = FingerprintRegistry()
```

**Implementation Time**: 3 days

---

### 1.5 Database Layer

#### SQLAlchemy Models

**File: `src/pubchem_finder/database/models.py`**
```python
"""Database models using SQLAlchemy."""
from datetime import datetime
from sqlalchemy import Column, Integer, String, Float, DateTime, Text, LargeBinary, JSON
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.dialects.postgresql import UUID
import uuid

Base = declarative_base()


class MoleculeModel(Base):
    """
    Relational storage for molecule metadata and properties.

    Fingerprints are stored separately in vector DB for efficient similarity search.
    This table stores metadata, properties, and references.
    """
    __tablename__ = "molecules"

    # Primary key
    id = Column(String(50), primary_key=True, default=lambda: str(uuid.uuid4()))

    # Chemical identifiers
    smiles = Column(String(500), nullable=False)
    canonical_smiles = Column(String(500), nullable=False, index=True)
    inchi = Column(Text, nullable=True)
    inchi_key = Column(String(27), nullable=True, index=True)  # Standard InChI key

    # Molecular properties (for filtering)
    molecular_weight = Column(Float, nullable=True)
    logp = Column(Float, nullable=True)
    num_h_donors = Column(Integer, nullable=True)
    num_h_acceptors = Column(Integer, nullable=True)
    tpsa = Column(Float, nullable=True)
    num_rotatable_bonds = Column(Integer, nullable=True)
    num_aromatic_rings = Column(Integer, nullable=True)

    # Metadata
    source = Column(String(50), default="pubchem", index=True)
    metadata = Column(JSON, default=dict)  # Flexible metadata storage

    # Timestamps
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow)

    # Indexes for common queries
    __table_args__ = (
        # Add composite indexes for property-based filtering
        # Index('ix_mol_weight_logp', 'molecular_weight', 'logp'),
    )


class FingerprintModel(Base):
    """
    Store fingerprints in relational DB (alternative to vector DB).

    For Phase 1, we'll store fingerprints here.
    Later, we can migrate to dedicated vector DB for better performance.
    """
    __tablename__ = "fingerprints"

    id = Column(Integer, primary_key=True, autoincrement=True)
    molecule_id = Column(String(50), nullable=False, index=True)
    fingerprint_type = Column(String(50), nullable=False, index=True)

    # Store as binary blob or array
    fingerprint_data = Column(LargeBinary, nullable=False)  # Numpy array serialized

    # Denormalized SMILES for quick reference
    smiles = Column(String(500), nullable=False)

    created_at = Column(DateTime, default=datetime.utcnow)

    __table_args__ = (
        # Unique constraint: one fingerprint type per molecule
        # UniqueConstraint('molecule_id', 'fingerprint_type'),
    )
```

#### Repository Pattern

**File: `src/pubchem_finder/database/repository.py`**
```python
"""Repository pattern for data access."""
from typing import Optional, List
from sqlalchemy.orm import Session
from sqlalchemy import and_
from .models import MoleculeModel, FingerprintModel
from ..core.molecule import Molecule
import numpy as np
import pickle


class MoleculeRepository:
    """
    Data access layer for molecules.

    Abstracts database operations from business logic.
    Makes it easy to swap database backends later.
    """

    def __init__(self, session: Session):
        self.session = session

    def save(self, molecule: Molecule) -> str:
        """Save or update molecule."""
        model = MoleculeModel(
            id=molecule.id or None,
            smiles=molecule.smiles,
            canonical_smiles=molecule.canonical_smiles,
            inchi=molecule.inchi,
            inchi_key=molecule.inchi_key,
            molecular_weight=molecule.molecular_weight,
            logp=molecule.logp,
            num_h_donors=molecule.num_h_donors,
            num_h_acceptors=molecule.num_h_acceptors,
            tpsa=molecule.tpsa,
            num_rotatable_bonds=molecule.num_rotatable_bonds,
            num_aromatic_rings=molecule.num_aromatic_rings,
            source=molecule.source,
            metadata=molecule.metadata,
        )

        self.session.add(model)
        self.session.flush()

        return model.id

    def save_fingerprint(self, molecule_id: str, fp_type: str,
                         fingerprint: np.ndarray, smiles: str):
        """Save fingerprint for a molecule."""
        fp_model = FingerprintModel(
            molecule_id=molecule_id,
            fingerprint_type=fp_type,
            fingerprint_data=pickle.dumps(fingerprint),  # Serialize numpy array
            smiles=smiles
        )
        self.session.add(fp_model)

    def get_by_id(self, molecule_id: str) -> Optional[Molecule]:
        """Retrieve molecule by ID."""
        model = self.session.query(MoleculeModel).filter(
            MoleculeModel.id == molecule_id
        ).first()

        if not model:
            return None

        return self._model_to_molecule(model)

    def get_by_smiles(self, smiles: str) -> Optional[Molecule]:
        """Retrieve molecule by canonical SMILES."""
        model = self.session.query(MoleculeModel).filter(
            MoleculeModel.canonical_smiles == smiles
        ).first()

        return self._model_to_molecule(model) if model else None

    def _model_to_molecule(self, model: MoleculeModel) -> Molecule:
        """Convert database model to domain object."""
        return Molecule(
            id=model.id,
            smiles=model.smiles,
            canonical_smiles=model.canonical_smiles,
            inchi=model.inchi,
            inchi_key=model.inchi_key,
            molecular_weight=model.molecular_weight,
            logp=model.logp,
            num_h_donors=model.num_h_donors,
            num_h_acceptors=model.num_h_acceptors,
            tpsa=model.tpsa,
            num_rotatable_bonds=model.num_rotatable_bonds,
            num_aromatic_rings=model.num_aromatic_rings,
            source=model.source,
            metadata=model.metadata,
            created_at=model.created_at,
            updated_at=model.updated_at,
        )

    def commit(self):
        """Commit transaction."""
        self.session.commit()

    def rollback(self):
        """Rollback transaction."""
        self.session.rollback()
```

**Implementation Time**: 3 days

---

### 1.6 PubChem Ingestion Pipeline

**File: `src/pubchem_finder/ingestion/pipeline.py`**
```python
"""Molecule ingestion pipeline."""
from pathlib import Path
from typing import Iterator
from rdkit import Chem
from ..chemistry.standardizer import MoleculeStandardizer
from ..chemistry.properties import calculate_properties
from ..core.molecule import Molecule
from ..fingerprints.registry import fingerprint_registry
from ..database.repository import MoleculeRepository
import logging

logger = logging.getLogger(__name__)


class IngestionPipeline:
    """
    Process molecules through standardization, fingerprinting, and storage.

    This is the main entry point for indexing molecules.
    """

    def __init__(
        self,
        repository: MoleculeRepository,
        fingerprint_types: list[str] = None,
    ):
        self.repository = repository
        self.standardizer = MoleculeStandardizer()
        self.fingerprint_types = fingerprint_types or ["morgan_2", "maccs"]

    def process_sdf(self, sdf_path: Path, batch_size: int = 1000):
        """
        Process SDF file and index molecules.

        Args:
            sdf_path: Path to SDF file
            batch_size: Number of molecules to process before commit
        """
        supplier = Chem.SDMolSupplier(str(sdf_path))

        processed = 0
        errors = 0

        for i, mol in enumerate(supplier):
            if mol is None:
                errors += 1
                continue

            try:
                # Process molecule
                molecule = self._process_molecule(mol)
                if molecule:
                    # Save to database
                    mol_id = self.repository.save(molecule)

                    # Generate and save fingerprints
                    self._save_fingerprints(mol_id, mol, molecule.canonical_smiles)

                    processed += 1

                    # Commit in batches
                    if processed % batch_size == 0:
                        self.repository.commit()
                        logger.info(f"Processed {processed} molecules")

            except Exception as e:
                logger.error(f"Error processing molecule {i}: {e}")
                errors += 1

        # Final commit
        self.repository.commit()
        logger.info(f"Completed: {processed} processed, {errors} errors")

    def _process_molecule(self, mol: Chem.Mol) -> Optional[Molecule]:
        """Process a single molecule."""
        # Standardize
        std_mol = self.standardizer.standardize(mol)
        if std_mol is None:
            return None

        # Get canonical SMILES
        canonical_smiles = Chem.MolToSmiles(std_mol)

        # Calculate properties
        props = calculate_properties(std_mol)

        # Create Molecule object
        molecule = Molecule(
            smiles=Chem.MolToSmiles(mol),  # Original
            canonical_smiles=canonical_smiles,
            inchi=Chem.MolToInchi(std_mol),
            inchi_key=Chem.MolToInchiKey(std_mol),
            **props
        )

        return molecule

    def _save_fingerprints(self, mol_id: str, mol: Chem.Mol, smiles: str):
        """Generate and save fingerprints."""
        for fp_type in self.fingerprint_types:
            generator = fingerprint_registry.get(fp_type)
            fingerprint = generator.generate(mol)
            self.repository.save_fingerprint(mol_id, fp_type, fingerprint, smiles)
```

**Implementation Time**: 2 days

---

### 1.7 Configuration Management

**File: `src/pubchem_finder/config/settings.py`**
```python
"""Application configuration using Pydantic Settings."""
from pydantic_settings import BaseSettings
from typing import List


class Settings(BaseSettings):
    """Application settings with environment variable support."""

    # Database
    database_url: str = "sqlite:///./pubchem.db"

    # Fingerprints
    default_fingerprints: List[str] = ["morgan_2", "maccs"]

    # Ingestion
    batch_size: int = 1000

    # Logging
    log_level: str = "INFO"

    class Config:
        env_file = ".env"
        env_file_encoding = "utf-8"


settings = Settings()
```

**Implementation Time**: 1 day

---

### 1.8 CLI Interface

**File: `src/pubchem_finder/cli/commands.py`**
```python
"""Command-line interface."""
import click
from pathlib import Path
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from ..config.settings import settings
from ..database.models import Base
from ..database.repository import MoleculeRepository
from ..ingestion.pipeline import IngestionPipeline


@click.group()
def cli():
    """PubChem Finder CLI."""
    pass


@cli.command()
def init_db():
    """Initialize database schema."""
    engine = create_engine(settings.database_url)
    Base.metadata.create_all(engine)
    click.echo("Database initialized successfully!")


@cli.command()
@click.argument("sdf_file", type=click.Path(exists=True))
@click.option("--batch-size", default=1000, help="Batch size for processing")
def index(sdf_file: str, batch_size: int):
    """Index molecules from SDF file."""
    engine = create_engine(settings.database_url)
    Session = sessionmaker(bind=engine)
    session = Session()

    repository = MoleculeRepository(session)
    pipeline = IngestionPipeline(repository)

    click.echo(f"Processing {sdf_file}...")
    pipeline.process_sdf(Path(sdf_file), batch_size=batch_size)
    click.echo("Indexing complete!")


@cli.command()
@click.argument("smiles")
def lookup(smiles: str):
    """Look up molecule by SMILES."""
    engine = create_engine(settings.database_url)
    Session = sessionmaker(bind=engine)
    session = Session()

    repository = MoleculeRepository(session)
    molecule = repository.get_by_smiles(smiles)

    if molecule:
        click.echo(f"Found: {molecule}")
    else:
        click.echo("Molecule not found")


if __name__ == "__main__":
    cli()
```

**Implementation Time**: 1 day

---

### 1.9 Testing Strategy

#### Unit Tests Structure
```
tests/
├── unit/
│   ├── test_standardizer.py
│   ├── test_fingerprints.py
│   ├── test_properties.py
│   └── test_repository.py
├── integration/
│   └── test_pipeline.py
└── fixtures/
    └── sample_molecules.sdf
```

**Example: `tests/unit/test_fingerprints.py`**
```python
"""Unit tests for fingerprint generators."""
import pytest
from rdkit import Chem
import numpy as np
from pubchem_finder.fingerprints.morgan import MorganFingerprint
from pubchem_finder.fingerprints.maccs import MACCSFingerprint


def test_morgan_fingerprint_shape():
    """Test Morgan fingerprint generates correct shape."""
    mol = Chem.MolFromSmiles("CCO")  # Ethanol
    generator = MorganFingerprint(radius=2, n_bits=2048)

    fp = generator.generate(mol)

    assert fp.shape == (2048,)
    assert fp.dtype == np.uint8


def test_morgan_fingerprint_deterministic():
    """Test fingerprint is deterministic."""
    mol = Chem.MolFromSmiles("CCO")
    generator = MorganFingerprint(radius=2, n_bits=2048)

    fp1 = generator.generate(mol)
    fp2 = generator.generate(mol)

    np.testing.assert_array_equal(fp1, fp2)


def test_maccs_fingerprint_size():
    """Test MACCS fingerprint is 167 bits."""
    mol = Chem.MolFromSmiles("c1ccccc1")  # Benzene
    generator = MACCSFingerprint()

    fp = generator.generate(mol)

    assert fp.shape == (167,)
```

**Implementation Time**: 3 days for comprehensive tests

---

### Phase 1 Deliverable Checklist

After Phase 1, you will have:

- [x] Clean, modular project structure
- [x] RDKit-based molecule processing
- [x] Morgan (ECFP) and MACCS fingerprints
- [x] Extensible fingerprint registry
- [x] SQLAlchemy database models
- [x] Repository pattern for data access
- [x] Molecule standardization pipeline
- [x] Property calculation
- [x] SDF ingestion pipeline
- [x] CLI for basic operations
- [x] Comprehensive unit tests (>80% coverage)
- [x] Documentation

**Total Implementation Time: 3-4 weeks**

---

## Phase 2: Basic Similarity Search (Weeks 5-7)

### Goal
Add similarity search capabilities with simple vector storage.

### Deliverables
✅ Similarity search implementation
✅ Multiple similarity metrics (Tanimoto, Dice)
✅ Basic vector storage (SQLite/PostgreSQL with pgvector)
✅ Search optimization (indexes)
✅ CLI search commands

### Key Components

#### 2.1 Similarity Search Engine
```python
# src/pubchem_finder/search/similarity.py
class SimilaritySearchEngine:
    """Simple brute-force similarity search."""

    def search(
        self,
        query_fingerprint: np.ndarray,
        fingerprint_type: str,
        metric: str = "tanimoto",
        threshold: float = 0.7,
        limit: int = 100,
    ) -> List[SearchResult]:
        """
        Search for similar molecules.

        For Phase 2, this will be brute-force.
        In Phase 3+, we'll optimize with vector DB.
        """
        pass
```

#### 2.2 Optional: pgvector Integration
```sql
-- Enable pgvector extension for PostgreSQL
CREATE EXTENSION vector;

-- Store fingerprints as vectors
ALTER TABLE fingerprints
ADD COLUMN fingerprint_vector vector(2048);

-- Create index for fast similarity search
CREATE INDEX ON fingerprints
USING ivfflat (fingerprint_vector vector_cosine_ops);
```

**Implementation Time: 2-3 weeks**

---

## Phase 3: REST API (Weeks 8-10)

### Goal
Build FastAPI-based REST API for programmatic access.

### Deliverables
✅ FastAPI application
✅ RESTful endpoints for search, lookup, indexing
✅ Request/response validation (Pydantic)
✅ API documentation (OpenAPI/Swagger)
✅ Rate limiting and authentication

### Key Endpoints
```python
# POST /api/v1/search/similarity
# GET  /api/v1/molecules/{id}
# POST /api/v1/molecules/batch
# GET  /api/v1/health
```

**Implementation Time: 2-3 weeks**

---

## Phase 4: Vector Database Migration (Weeks 11-13)

### Goal
Migrate to dedicated vector database for performance.

### Deliverables
✅ Qdrant/Milvus integration
✅ Data migration tools
✅ Performance benchmarks
✅ Hybrid search (properties + similarity)

**Implementation Time: 2-3 weeks**

---

## Future Phases (Post-Foundation)

### Phase 5: Machine Learning Integration
- GNN-based embeddings
- Pre-trained molecular models
- Hybrid similarity scoring

### Phase 6: Advanced Features
- Substructure search
- Scaffold hopping
- Batch operations
- Real-time search

### Phase 7: Production Hardening
- Kubernetes deployment
- Monitoring and observability
- Load balancing
- Caching layer

---

## Technology Decisions

### Phase 1 Tech Stack
```toml
[tool.poetry.dependencies]
python = "^3.10"
rdkit = "^2023.9.1"
sqlalchemy = "^2.0.0"
alembic = "^1.13.0"
pydantic = "^2.5.0"
pydantic-settings = "^2.1.0"
click = "^8.1.7"
loguru = "^0.7.2"

[tool.poetry.group.dev.dependencies]
pytest = "^7.4.0"
pytest-cov = "^4.1.0"
black = "^23.11.0"
ruff = "^0.1.6"
mypy = "^1.7.0"
```

### Database Choice for Phase 1
**SQLite** for development, **PostgreSQL** for production.

Why?
- SQLite: Zero setup, perfect for development and testing
- PostgreSQL: Production-ready, supports pgvector for Phase 2+

### Future: Vector DB Choice
**Qdrant** (recommended) or **Milvus**

Why Qdrant?
- Easier setup than Milvus
- Excellent Python SDK
- Good performance for <10M molecules
- Hybrid search support

---

## Development Workflow

### Branch Strategy
```
main (production-ready code)
├── develop (integration branch)
│   ├── feature/phase1-core
│   ├── feature/phase1-fingerprints
│   ├── feature/phase1-ingestion
│   └── feature/phase2-search
```

### CI/CD Pipeline
```yaml
# .github/workflows/ci.yml
name: CI
on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      - uses: actions/setup-python@v4
        with:
          python-version: '3.10'
      - run: pip install poetry
      - run: poetry install
      - run: poetry run pytest --cov
      - run: poetry run black --check .
      - run: poetry run ruff check .
      - run: poetry run mypy .
```

---

## Success Metrics - Phase 1

### Technical Metrics
- [ ] Unit test coverage >80%
- [ ] Type hint coverage >90%
- [ ] Can process 10,000 molecules/minute
- [ ] Zero crashes on valid input

### Functional Metrics
- [ ] Successfully indexes PubChem sample dataset
- [ ] Fingerprints match RDKit reference implementation
- [ ] Standardization produces consistent SMILES

### Code Quality
- [ ] All code passes linting (ruff)
- [ ] All code formatted (black)
- [ ] Type checking passes (mypy)
- [ ] Documentation coverage >70%

---

## Risk Mitigation

| Risk | Mitigation |
|------|------------|
| RDKit installation issues | Provide Docker dev environment |
| Large SDF files crash | Implement streaming parser |
| Database schema changes | Use Alembic migrations |
| Fingerprint bugs | Extensive unit tests vs. RDKit reference |

---

## Getting Started

### Day 1-2: Project Setup
```bash
# Create new branch
git checkout -b feature/phase1-foundation

# Initialize project with Poetry
poetry init
poetry add rdkit sqlalchemy pydantic click
poetry add --group dev pytest black ruff mypy

# Create directory structure
mkdir -p src/pubchem_finder/{core,chemistry,fingerprints,database,ingestion,config,cli}
mkdir -p tests/{unit,integration,fixtures}

# Initialize git hooks
poetry add --group dev pre-commit
pre-commit install
```

### Week 1: Core + Chemistry
Focus on `core/` and `chemistry/` modules with tests.

### Week 2: Fingerprints
Implement fingerprint generators with comprehensive tests.

### Week 3: Database + Ingestion
Build database layer and ingestion pipeline.

### Week 4: Integration + Testing
End-to-end testing, documentation, bug fixes.

---

## Conclusion

This plan prioritizes:

1. **Solid Foundation**: Well-tested, clean code
2. **Extensibility**: Plugin architecture for future features
3. **Pragmatism**: Start simple (SQLite), scale later (PostgreSQL + vector DB)
4. **Incremental Value**: Each phase delivers working software

By end of Phase 1, you'll have a production-quality molecular indexing system that can be extended with ML, advanced search, and APIs without rewriting core components.

**Next Step**: Start with Phase 1, Week 1 - Core & Chemistry modules!
