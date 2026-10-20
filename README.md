# PubChem Finder

Molecular similarity search over [PubChem](https://pubchem.ncbi.nlm.nih.gov/), built on RDKit.
The goal is a tool you can run on a laptop: build a local fingerprint index of PubChem
Compound, search it from the command line, and expose it to Claude (or any MCP client)
as an [MCP](https://modelcontextprotocol.io/) server.

[![CI](https://github.com/MysterionRise/pubchem-finder/actions/workflows/ci.yml/badge.svg)](https://github.com/MysterionRise/pubchem-finder/actions/workflows/ci.yml)
[![Python 3.11+](https://img.shields.io/badge/python-3.11%2B-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

> **Status:** being rebuilt. The chemistry core (standardization, descriptors,
> fingerprints, similarity metrics) is in place and tested. Index building, search,
> the CLI and the MCP server are the next milestones. See [ROADMAP.md](ROADMAP.md).

## What is here today

- `pubchem_finder.chemistry`: `MoleculeStandardizer` (RDKit `rdMolStandardize`:
  cleanup, largest organic fragment, neutralize) and `calculate_properties`
  (MW, LogP, HBD/HBA, TPSA, rotatable bonds, aromatic rings) plus `is_druglike`.
- `pubchem_finder.fingerprints`: `MorganFingerprint` (ECFP/FCFP via
  `rdFingerprintGenerator`), `MACCSFingerprint`, `TanimotoSimilarity`,
  `DiceSimilarity`, `CosineSimilarity`, `bulk_tanimoto`, and a small registry
  (`ecfp4`, `ecfp6`, `fcfp4`, `maccs`, ...).
- `pubchem_finder.core`: the `Molecule` record and the abstract
  `FingerprintGenerator` / `SimilarityMetric` interfaces.

```python
from pubchem_finder.chemistry import MoleculeStandardizer, calculate_properties
from pubchem_finder.fingerprints import MorganFingerprint, TanimotoSimilarity

std = MoleculeStandardizer()
print(std.standardize_smiles("CC(=O)[O-].[Na+]"))  # CC(=O)O

gen = MorganFingerprint(radius=2, n_bits=1024)  # ECFP4, 1024 bits
aspirin = gen.generate_from_smiles("CC(=O)Oc1ccccc1C(=O)O")
salicylic = gen.generate_from_smiles("O=C(O)c1ccccc1O")
print(TanimotoSimilarity().calculate(aspirin, salicylic))
```

## Development

Requires [uv](https://docs.astral.sh/uv/). Python 3.11, 3.12 and 3.13 are supported.

```bash
uv sync            # create .venv and install everything, including dev tools
make check         # ruff lint + format check, mypy (strict), pytest with coverage
make test          # just the tests
```

Pre-commit hooks (ruff lint and format) can be installed with `uv run pre-commit install`
after `uv tool install pre-commit`.

## Planned architecture

```
PubChem Compound SDF (HTTPS) --> streaming RDKit parse --> Parquet metadata (CID, SMILES, InChIKey, MW, formula)
                                                      \-> FPSim2 fingerprint index (ECFP4, 1024 bits)
pubchem-finder CLI  ----\
                         +--> search facade --> FPSim2 (exact Tanimoto, in-memory or on-disk) + DuckDB
MCP server (stdio) -----/                   \-> PubChem PUG REST for live lookups without an index
```

The index engine is [FPSim2](https://github.com/chembl/FPSim2), the similarity
search library behind ChEMBL's interface. A PostgreSQL + RDKit cartridge backend for
substructure search is an optional later addition behind the same interface.

## History

This repository started in 2020 as an Indigo + Elasticsearch indexer for PubChem
dumps. That code has been removed; it is available in the git history before
commit `06619d9` if you need it.

## License

MIT. See [LICENSE](LICENSE).
