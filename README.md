# PubChem Finder

Molecular similarity search over [PubChem](https://pubchem.ncbi.nlm.nih.gov/), built on RDKit.
The goal is a tool you can run on a laptop: build a local fingerprint index of PubChem
Compound, search it from the command line, and expose it to Claude (or any MCP client)
as an [MCP](https://modelcontextprotocol.io/) server.

[![CI](https://github.com/MysterionRise/pubchem-finder/actions/workflows/ci.yml/badge.svg)](https://github.com/MysterionRise/pubchem-finder/actions/workflows/ci.yml)
[![Python 3.11+](https://img.shields.io/badge/python-3.11%2B-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

> **Status:** PubChem lookups and live structure search, synonyms, cross-references,
> curated safety and pharmacology annotations, bioassay summaries, RDKit descriptors,
> standardization, fingerprint comparison and a local similarity index over PubChem
> (or any SDF files) all work from the CLI and from an MCP server.
> See [ROADMAP.md](ROADMAP.md) for what comes next.

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

### Similarity search over a local index

Build an index once, then search it as often as you like, offline and without rate limits.
The examples assume `pubchem-finder` is installed with `uv tool install` as above.

```bash
pubchem-finder index build --limit 1        # first PubChem chunk: ~500k compounds, a few minutes
pubchem-finder search "CC(=O)Oc1ccccc1C(=O)O" --threshold 0.7 --top 20
pubchem-finder index status

pubchem-finder index build --sdf my_library.sdf.gz   # or index your own SDF files
```

### Ask PubChem directly

The `pubchem` commands query PubChem's live services, with no index needed:

```bash
pubchem-finder pubchem similarity "CC(=O)Oc1ccccc1C(=O)O" --threshold 90   # across all ~120M compounds
pubchem-finder pubchem substructure "c1ccc2[nH]ccc2c1" --max 50              # everything containing indole
pubchem-finder pubchem exact "CC(=O)Oc1ccccc1C(=O)O" --identity same_connectivity
pubchem-finder pubchem synonyms 2244          # names and CAS numbers
pubchem-finder pubchem xrefs 2244             # registry IDs, PubMed articles, patents
pubchem-finder pubchem describe 2244          # curated descriptions with sources
pubchem-finder pubchem annotations 2244 --heading "GHS Classification"   # safety, pharmacology, ...
pubchem-finder pubchem bioactivity 2244       # assay outcomes, top targets, most potent results
```

PubChem's similarity search uses its own fingerprint and does not report scores, so each hit is
also scored here with ECFP4, the same fingerprint as the local index, and sorted by it. Annotations
name the organisation that deposited them; some carry their own license terms.

Successful PubChem responses are cached for 7 days in the data directory, so repeated questions
cost nothing and stay within PubChem's rate limits. `pubchem-finder cache info` and
`pubchem-finder cache clear` manage it; set `PUBCHEM_FINDER_CACHE=0` to turn it off.

### Pictures, triage and batch files

```bash
# Draw structures; highlight a SMARTS pattern, what each analogue adds to a reference,
# or which atoms drive ECFP4 similarity to it (similarity map)
pubchem-finder draw "CC(=O)Oc1ccccc1C(=O)O" "COC(=O)c1ccccc1O" -l aspirin -l "methyl salicylate" \
    --highlight diff --reference "CC(=O)Oc1ccccc1C(=O)O" -o analogues.png

# Med-chem triage: QED, synthetic accessibility, Lipinski/Veber, PAINS/Brenk/NIH alerts, scaffold
pubchem-finder triage "O=C1C(=Cc2ccccc2)SC(=S)N1"

# Whole files: .smi ("SMILES [id]" per line), .csv/.tsv with a smiles column, or - for stdin
pubchem-finder triage -i library.smi -o triage.csv
pubchem-finder props -i library.csv -o props.csv
pubchem-finder standardize -i library.smi -o clean.csv
pubchem-finder search -i queries.smi -o hits.csv --threshold 0.6 --top 50   # one row per hit
pubchem-finder search "CC(=O)Oc1ccccc1C(=O)O" --group-by-scaffold
```

In `diff` mode each drawing is aligned to the core it shares with the reference, and only the
atoms the analogue adds are highlighted; groups the analogue lacks cannot be drawn on it. Rows
that fail in batch mode get an `error` column instead of stopping the run.

`index build --limit N` downloads the first N of PubChem's ~350 SDF chunks (about 500,000
CIDs each) from the NCBI HTTPS mirror, verifies each file's MD5, fingerprints it, deletes the
download, and merges everything into one searchable index. Re-running is safe: finished chunks
are skipped, so an interrupted build resumes and a larger `--limit` extends an existing index.
`--all` indexes all of PubChem. The default fingerprint is ECFP4 (Morgan radius 2, 1024 bits);
`--radius` and `--fp-size` change it for a new index.

Data lives in the platform data directory (`~/.local/share/pubchem-finder` on Linux,
`~/Library/Application Support/pubchem-finder` on macOS, `%LOCALAPPDATA%\pubchem-finder` on
Windows). Override it with `--home DIR` or `PUBCHEM_FINDER_HOME`.

**How big is it?** The first column is a real run of `index build --limit 1` on a MacBook Pro.
The full-PubChem column scales those numbers to PubChem's roughly 120 million compounds and
is an estimate, not a measurement.

| | First PubChem chunk (measured) | All of PubChem (estimated) |
|---|---|---|
| Compounds indexed | 444,794 (29 skipped by RDKit) | ~120 million |
| Fingerprinting | 80 s on one core (~5,500 compounds/s) | ~1 h on 8 cores, plus download time |
| Download | one SDF file, deleted after indexing | ~100 GB over HTTPS, streamed chunk by chunk |
| Index file on disk | 36 MB | ~10 GB |
| RAM for in-memory search | 64 MB | ~17 GB |
| Free disk while building | small | ~45 GB (chunk files, merged index, sort scratch, metadata) |

Search latency, measured on 100,000 compounds in a 4-core Linux sandbox: about 2 ms per
query in memory and 210 ms on disk. Both grow roughly linearly with index size, so a
full-PubChem query takes about 1-3 s in memory and minutes on disk.

The index is searched in memory when it fits comfortably in available RAM and from disk
otherwise; set `PUBCHEM_FINDER_INDEX_MODE=memory|disk` to force a mode. On a machine without
~20 GB of free RAM, index a subset with `--limit` rather than all of PubChem.

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
| `similarity_search(smiles, threshold, top_k, metric)` | Most similar compounds in the local index, with CID, SMILES, IUPAC name, formula, MW and PubChem link. Needs `index build`. |
| `index_status()` | Whether the local index is built, its size and fingerprint, and how to build it. |
| `pubchem_similarity_search(smiles, threshold, max_records)` | Similar compounds across all of PubChem, live, scored with ECFP4 and sorted. |
| `pubchem_substructure_search(smiles, max_records)` | Compounds across all of PubChem that contain a fragment. |
| `pubchem_exact_match(smiles, identity_type)` | The PubChem record(s) for exactly this structure, optionally ignoring stereo or isotopes. |
| `pubchem_synonyms(cid, limit)` | Common, brand and systematic names, plus CAS Registry Numbers. |
| `pubchem_xrefs(cid)` | Registry IDs (CAS, ChEMBL, DrugBank...), PubMed articles and patents. |
| `pubchem_description(cid)` | Curated descriptions with sources. |
| `pubchem_annotations(cid, heading, max_chars)` | One section of the compound page: GHS classification, hazards, pharmacology, mechanism of action, drug information, experimental properties, toxicity, use, food additives. Every item names its source. |
| `pubchem_bioactivity(cid)` | Bioassay summary: assay count, outcomes, top targets, most potent active results. |
| `draw_molecules(smiles, legends, highlight, smarts, reference_smiles)` | Structures as an inline image, with SMARTS, "what changed" or similarity-map highlighting. Offline. |
| `medchem_profile(smiles)` | QED, synthetic accessibility, Lipinski/Veber, PAINS/Brenk/NIH alerts with atoms, Murcko scaffold. Offline. |

`similarity_search` also takes `group_by_scaffold` to show how many chemical series the hits span.

Then ask Claude things like "what is the molecular weight and XLogP of caffeine?",
"how similar are aspirin and salicylic acid?", "is ibuprofen hazardous, and what does its
GHS label say?", "what are the main targets of celecoxib in PubChem assays?", "show me the
closest analogues of aspirin and highlight what changed", "does this compound have any PAINS
alerts?" or, with an index built, "find the 10 compounds most similar to ibuprofen".

The server loads the index on the first search, not at startup, so it starts instantly. Point it
at a non-default data directory by adding `"env": {"PUBCHEM_FINDER_HOME": "/path"}` to the
Claude Desktop entry, or `--env PUBCHEM_FINDER_HOME=/path` to `claude mcp add`.

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

Requires [uv](https://docs.astral.sh/uv/). Python 3.11, 3.12 and 3.13 are supported; the
repository's `.python-version` makes `uv sync` use 3.12 even if pyenv or conda put an older
Python first on your PATH.

```bash
uv sync            # create .venv and install everything, including dev tools
make check         # ruff lint + format check, mypy (strict), pytest with coverage
make test          # just the tests
PUBCHEM_LIVE=1 uv run pytest -m live   # tests that call the real PubChem services
```

Pre-commit hooks (ruff lint and format) can be installed with `uv run pre-commit install`
after `uv tool install pre-commit`.

## Architecture

```
PubChem Compound SDF (HTTPS) --> streaming RDKit parse --> chunks/<chunk>.parquet  (CID, SMILES, InChIKey, name, formula, MW)
                                                      \-> chunks/<chunk>.h5       (FPSim2 fingerprints)
                                                            \-> merge + popcount sort -> index/pubchem.h5 + manifest.json
pubchem-finder CLI  ----\
                         +--> service --> FPSim2 (exact similarity, in memory or on disk) + DuckDB over Parquet
MCP server (stdio) -----/            \-> PubChem PUG REST for live lookups
```

The index engine is [FPSim2](https://github.com/chembl/FPSim2), the similarity search
library behind ChEMBL's web interface. Each SDF chunk is fingerprinted independently, so
builds run in parallel and resume after interruption.

## History

This repository started in 2020 as an Indigo + Elasticsearch indexer for PubChem
dumps. That code has been removed; it is available in the git history before
commit `06619d9` if you need it.

## License

MIT. See [LICENSE](LICENSE).
