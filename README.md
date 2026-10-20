# PubChem Finder

Molecular similarity search over [PubChem](https://pubchem.ncbi.nlm.nih.gov/), built on RDKit.
The goal is a tool you can run on a laptop: build a local fingerprint index of PubChem
Compound, search it from the command line, and expose it to Claude (or any MCP client)
as an [MCP](https://modelcontextprotocol.io/) server.

[![CI](https://github.com/MysterionRise/pubchem-finder/actions/workflows/ci.yml/badge.svg)](https://github.com/MysterionRise/pubchem-finder/actions/workflows/ci.yml)
[![Python 3.11+](https://img.shields.io/badge/python-3.11%2B-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

> **Status:** v0.1 of the rebuild. PubChem lookups, RDKit descriptors,
> standardization and fingerprint comparison work from the CLI and from an MCP
> server. The local similarity index over all of PubChem is the next milestone.
> See [ROADMAP.md](ROADMAP.md).

## Quick start

Requires [uv](https://docs.astral.sh/uv/). No install step: `uvx` fetches and caches the tool.

```bash
# Look a compound up in PubChem (by name, CID, SMILES or InChIKey)
uvx --from git+https://github.com/MysterionRise/pubchem-finder pubchem-finder lookup caffeine
uvx --from git+https://github.com/MysterionRise/pubchem-finder pubchem-finder lookup 2244 --by cid --json

# RDKit descriptors, standardization, similarity between two structures (offline)
uvx --from git+https://github.com/MysterionRise/pubchem-finder pubchem-finder props "CC(=O)Oc1ccccc1C(=O)O"
uvx --from git+https://github.com/MysterionRise/pubchem-finder pubchem-finder standardize "CC(=O)[O-].[Na+]"
uvx --from git+https://github.com/MysterionRise/pubchem-finder pubchem-finder compare "CCO" "CCCO"
```

Or install once: `uv tool install git+https://github.com/MysterionRise/pubchem-finder`, then use `pubchem-finder ...`.

### Use it from Claude (MCP)

The same operations are available as an [MCP](https://modelcontextprotocol.io/) server over stdio.

Claude Code:

```bash
claude mcp add pubchem-finder -- uvx --from git+https://github.com/MysterionRise/pubchem-finder pubchem-finder mcp
```

Claude Desktop (`claude_desktop_config.json`):

```json
{
  "mcpServers": {
    "pubchem-finder": {
      "command": "uvx",
      "args": ["--from", "git+https://github.com/MysterionRise/pubchem-finder", "pubchem-finder", "mcp"]
    }
  }
}
```

Tools exposed:

| Tool | What it does |
|---|---|
| `pubchem_lookup(query, by, limit)` | Resolve a name, CID, SMILES or InChIKey to PubChem records (identifiers, formula, MW, XLogP, H-bond counts, TPSA, link). Live PubChem PUG REST. |
| `compound_properties(smiles)` | RDKit descriptors and Lipinski check for any SMILES. Offline. |
| `standardize(smiles)` | Remove salts and counter-ions, neutralize, canonicalize. Offline. |
| `compare(smiles_a, smiles_b)` | Tanimoto similarity on ECFP4 and MACCS fingerprints. Offline. |

Then ask Claude things like "what is the molecular weight and XLogP of caffeine?" or
"how similar are aspirin and salicylic acid?".

## Library

The CLI and MCP server are thin layers over `pubchem_finder.service`; the chemistry
lives in `pubchem_finder.chemistry` and `pubchem_finder.fingerprints`.

```python
from pubchem_finder import service

service.properties("CC(=O)Oc1ccccc1C(=O)O").molecular_weight  # 180.16
service.standardize("CC(=O)[O-].[Na+]").smiles  # 'CC(=O)O'
service.compare("CCO", "CCCO").ecfp4_tanimoto  # 0.556
service.lookup("caffeine").compounds[0].cid  # 2519 (network)
```

Lower level:

```python
from pubchem_finder.chemistry import MoleculeStandardizer, calculate_properties
from pubchem_finder.fingerprints import MorganFingerprint, TanimotoSimilarity

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
PUBCHEM_LIVE=1 uv run pytest -m live   # the one test that calls the real PubChem API
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
