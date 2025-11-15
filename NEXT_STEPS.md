# Next Steps - Getting Started with Phase 1

## What We've Accomplished

✅ **Comprehensive Codebase Analysis** - Analyzed current Indigo/Elasticsearch implementation
✅ **Replatforming Strategy** - Created detailed technical specification in `REPLATFORMING_PROMPT.md`
✅ **Implementation Plan** - Designed 4-phase approach in `IMPLEMENTATION_PLAN.md`

---

## Quick Start Guide

### Start with Phase 1: Core Foundation (3-4 weeks)

Phase 1 builds the essential components that everything else depends on. No APIs, no ML, no complex features - just solid, well-tested foundation.

---

## Week 1: Project Setup & Core (Days 1-7)

### Day 1-2: Initialize Project

```bash
# Create feature branch
git checkout -b feature/phase1-foundation

# Initialize with Poetry (recommended) or PDM
poetry init
poetry add rdkit sqlalchemy alembic pydantic pydantic-settings click loguru
poetry add --group dev pytest pytest-cov black ruff mypy pre-commit

# Create directory structure
mkdir -p src/pubchem_finder/{core,chemistry,fingerprints,database,ingestion,config,cli}
mkdir -p tests/{unit,integration,fixtures}
touch src/pubchem_finder/__init__.py

# Set up pre-commit hooks
cat > .pre-commit-config.yaml <<'EOF'
repos:
  - repo: https://github.com/psf/black
    rev: 23.11.0
    hooks:
      - id: black
  - repo: https://github.com/astral-sh/ruff-pre-commit
    rev: v0.1.6
    hooks:
      - id: ruff
EOF

pre-commit install
```

### Day 3-4: Core Data Models

Implement:
- `src/pubchem_finder/core/molecule.py` - Molecule dataclass
- `src/pubchem_finder/core/fingerprint.py` - Abstract base classes
- `src/pubchem_finder/core/exceptions.py` - Custom exceptions

See `IMPLEMENTATION_PLAN.md` Section 1.2 for complete code.

**Tests**: `tests/unit/test_molecule.py`

### Day 5-7: Chemistry Layer

Implement:
- `src/pubchem_finder/chemistry/standardizer.py` - Molecule standardization
- `src/pubchem_finder/chemistry/properties.py` - Property calculation
- `src/pubchem_finder/chemistry/validator.py` - Validation

See `IMPLEMENTATION_PLAN.md` Section 1.3 for complete code.

**Tests**: `tests/unit/test_standardizer.py`, `tests/unit/test_properties.py`

---

## Week 2: Fingerprints (Days 8-14)

### Implement Fingerprint Generators

Priority order:
1. **Morgan/ECFP** (most important) - `src/pubchem_finder/fingerprints/morgan.py`
2. **MACCS keys** - `src/pubchem_finder/fingerprints/maccs.py`
3. **RDKit topological** (optional) - `src/pubchem_finder/fingerprints/rdkit_fp.py`

### Implement Registry Pattern

- `src/pubchem_finder/fingerprints/registry.py` - Central fingerprint registry

This makes adding new fingerprints (including ML-based ones later) trivial!

### Implement Similarity Metrics

- `src/pubchem_finder/fingerprints/base.py` - Tanimoto, Dice, Cosine

**Tests**: Extensive tests in `tests/unit/test_fingerprints.py`

See `IMPLEMENTATION_PLAN.md` Section 1.4 for complete code.

---

## Week 3: Database & Ingestion (Days 15-21)

### Database Layer

1. Define SQLAlchemy models - `src/pubchem_finder/database/models.py`
2. Implement repository pattern - `src/pubchem_finder/database/repository.py`
3. Set up Alembic migrations

```bash
poetry run alembic init alembic
# Edit alembic.ini and alembic/env.py
poetry run alembic revision --autogenerate -m "Initial schema"
```

### Ingestion Pipeline

Implement `src/pubchem_finder/ingestion/pipeline.py`:
- SDF file parsing
- Molecule standardization
- Fingerprint generation
- Batch database insertion

**Tests**: Integration tests in `tests/integration/test_pipeline.py`

See `IMPLEMENTATION_PLAN.md` Sections 1.5-1.6 for complete code.

---

## Week 4: CLI & Polish (Days 22-28)

### CLI Commands

Implement `src/pubchem_finder/cli/commands.py`:

```bash
pubchem-finder init-db           # Initialize database
pubchem-finder index <sdf_file>  # Index molecules
pubchem-finder lookup <smiles>   # Lookup by SMILES
pubchem-finder stats             # Show statistics
```

### Configuration

- `src/pubchem_finder/config/settings.py` - Pydantic settings with .env support

### Testing & Documentation

- Achieve 80%+ test coverage
- Write docstrings for all public APIs
- Create README with usage examples
- Fix any bugs found during integration testing

---

## Phase 1 Acceptance Criteria

Before moving to Phase 2, ensure:

- [ ] All unit tests pass with >80% coverage
- [ ] Successfully index a sample PubChem dataset (e.g., 10,000 molecules)
- [ ] Fingerprints validated against RDKit reference
- [ ] No memory leaks during large ingestion
- [ ] Code passes all linters (black, ruff, mypy)
- [ ] Documentation complete

---

## Technology Choices for Phase 1

### Packaging: Poetry vs PDM

**Recommended: Poetry**
- More mature, better IDE support
- Widely used in data science community

```bash
# Install Poetry
curl -sSL https://install.python-poetry.org | python3 -
```

### Database: SQLite vs PostgreSQL

**Start with SQLite**, migrate to PostgreSQL later:

```python
# Development
DATABASE_URL=sqlite:///./pubchem_dev.db

# Production (Phase 2+)
DATABASE_URL=postgresql://user:pass@localhost/pubchem
```

### Python Version

**Minimum: Python 3.10** (for better type hints)
**Recommended: Python 3.11+** (performance improvements)

---

## Quick Reference: File Checklist

Core files to create in Phase 1:

```
src/pubchem_finder/
├── core/
│   ├── __init__.py
│   ├── molecule.py          ✅ Week 1
│   ├── fingerprint.py       ✅ Week 1
│   └── exceptions.py        ✅ Week 1
│
├── chemistry/
│   ├── __init__.py
│   ├── standardizer.py      ✅ Week 1
│   ├── properties.py        ✅ Week 1
│   └── validator.py         ✅ Week 1
│
├── fingerprints/
│   ├── __init__.py
│   ├── base.py              ✅ Week 2
│   ├── morgan.py            ✅ Week 2
│   ├── maccs.py             ✅ Week 2
│   └── registry.py          ✅ Week 2
│
├── database/
│   ├── __init__.py
│   ├── models.py            ✅ Week 3
│   └── repository.py        ✅ Week 3
│
├── ingestion/
│   ├── __init__.py
│   ├── pipeline.py          ✅ Week 3
│   └── sdf_parser.py        ✅ Week 3
│
├── config/
│   ├── __init__.py
│   └── settings.py          ✅ Week 4
│
└── cli/
    ├── __init__.py
    └── commands.py          ✅ Week 4
```

---

## After Phase 1: What's Next?

### Phase 2: Basic Similarity Search (Weeks 5-7)

Once foundation is solid:
1. Implement similarity search engine
2. Add PostgreSQL with pgvector extension
3. Build CLI search commands
4. Benchmark performance

### Phase 3: REST API (Weeks 8-10)

Build FastAPI application:
- Similarity search endpoint
- Molecule lookup
- Batch operations
- OpenAPI docs

### Phase 4: Vector Database (Weeks 11-13)

Production-scale search:
- Migrate to Qdrant or Milvus
- Optimize for millions of molecules
- Hybrid search capabilities

### Future: ML Integration

Once core is working:
- Add GNN-based embeddings
- Integrate pre-trained models (ChemBERTa)
- Deep metric learning for similarity

---

## Resources

### Documentation
- [RDKit Getting Started](https://www.rdkit.org/docs/GettingStartedInPython.html)
- [SQLAlchemy Tutorial](https://docs.sqlalchemy.org/en/20/tutorial/)
- [Pydantic Documentation](https://docs.pydantic.dev/)
- [Poetry Documentation](https://python-poetry.org/docs/)

### Example Datasets
- [PubChem Sample Data](https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/SDF/)
- Use small files (e.g., first 10k compounds) for testing

### Code Examples
See `IMPLEMENTATION_PLAN.md` for:
- Complete code for all modules
- Testing strategies
- Configuration examples

---

## Development Tips

### 1. Test-Driven Development

Write tests FIRST for each module:
```python
# tests/unit/test_standardizer.py
def test_standardize_removes_salts():
    """Test that standardization removes salt counterions."""
    smiles = "CCO.Cl"  # Ethanol with chloride salt
    standardizer = MoleculeStandardizer()
    result = standardizer.standardize_smiles(smiles)
    assert result == "CCO"  # Just ethanol
```

Then implement the feature.

### 2. Use Type Hints Everywhere

```python
def standardize(self, mol: Chem.Mol) -> Optional[Chem.Mol]:
    """Type hints help catch bugs early."""
    pass
```

Run `mypy` frequently to catch type errors.

### 3. Commit Often

```bash
git add src/pubchem_finder/core/
git commit -m "Add core molecule data structures"

git add tests/unit/test_molecule.py
git commit -m "Add molecule tests"
```

Small, focused commits are easier to review and debug.

### 4. Use Docker for Development

Create `docker-compose.yml` for consistent environment:
```yaml
version: '3.8'
services:
  dev:
    image: continuumio/miniconda3
    volumes:
      - .:/app
    working_dir: /app
    command: bash
```

### 5. Benchmark Early

Create `scripts/benchmark.py` to track performance:
- Molecules processed per second
- Memory usage
- Fingerprint generation speed

---

## Questions to Consider

### Before Starting

1. **Target Scale**: How many molecules do you plan to index?
   - <1M: SQLite + simple search OK
   - 1M-10M: PostgreSQL + pgvector
   - 10M+: Dedicated vector DB from Phase 2

2. **Performance Requirements**: Query latency SLA?
   - <1s: Simple similarity search OK
   - <100ms: Need vector DB + optimization

3. **Use Cases**: What searches are most important?
   - Similarity search only: Focus on fingerprints
   - Substructure search needed: Add in Phase 3
   - ML-based ranking: Plan for Phase 5

---

## Support & Community

### RDKit Community
- [RDKit Mailing List](https://sourceforge.net/p/rdkit/mailman/)
- [GitHub Discussions](https://github.com/rdkit/rdkit/discussions)

### Cheminformatics Resources
- [Open-Source Cheminformatics Blog](https://www.cheminformania.com/)
- [Practical Cheminformatics](https://practicalcheminformatics.blogspot.com/)

---

## Let's Start Building! 🚀

**Immediate Action**: Begin with Week 1, Day 1-2 project setup.

```bash
# Right now:
git checkout -b feature/phase1-foundation
poetry init

# Start coding!
```

Good luck with the implementation! Remember: **Foundation first, features later.**
