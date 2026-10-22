# PubChem Finder

Search [PubChem](https://pubchem.ncbi.nlm.nih.gov/) by chemical structure from your
terminal or from Claude. PubChem Finder builds a local fingerprint index of PubChem (all
of it, or the slice you need) for exact similarity and substructure search with no rate
limits, wraps PubChem's live services for names, safety data and bioassays, and adds RDKit
chemistry: structure drawings, med-chem triage and descriptors. It runs as a command-line
tool, a Python library and an [MCP](https://modelcontextprotocol.io/) server for Claude
Desktop, Claude Code and other MCP clients.

[![PyPI](https://img.shields.io/pypi/v/pubchem-finder.svg)](https://pypi.org/project/pubchem-finder/)
[![Python](https://img.shields.io/pypi/pyversions/pubchem-finder.svg)](https://pypi.org/project/pubchem-finder/)
[![CI](https://github.com/MysterionRise/pubchem-finder/actions/workflows/ci.yml/badge.svg)](https://github.com/MysterionRise/pubchem-finder/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/MysterionRise/pubchem-finder/blob/master/LICENSE)

![A similarity search for aspirin and a SMARTS substructure search over 444,794 PubChem compounds, with the hits appearing and their key atoms highlighted](https://raw.githubusercontent.com/MysterionRise/pubchem-finder/master/docs/media/search-demo.gif)

*Real searches over the first 444,794 compounds of PubChem, rendered from the index the
project's real-data workflow builds: similarity hits for aspirin with the atoms that differ
from it highlighted, then hits for an aryl sulfonamide SMARTS with the match highlighted.
[MP4 version](https://raw.githubusercontent.com/MysterionRise/pubchem-finder/master/docs/media/search-demo.mp4).*

## Why this one

Other PubChem tools, MCP servers included, forward questions to PubChem's web API. That
works for looking compounds up, but structure search through the API is rate-limited,
times out on generic queries and returns unscored samples. PubChem Finder also keeps a
local index:

- **Exact, scored similarity search** over every indexed compound (ECFP4 Tanimoto by
  default) in milliseconds, using [FPSim2](https://github.com/chembl/FPSim2), the engine
  behind ChEMBL's similarity search.
- **Local substructure search** with SMILES or SMARTS: a fingerprint screen, then an
  atom-by-atom match, so every hit is a real match.
- **Offline chemistry** with RDKit: descriptors, standardization, PAINS/Brenk/NIH alerts,
  QED, synthetic accessibility, and drawings that show what changed between analogues.
- **PubChem knowledge, cited:** synonyms and CAS numbers, cross-references, GHS
  classification, pharmacology, toxicity and bioassay summaries, each with its source.

No database server, no Docker: the index is a pair of files in your data directory.

## Use it from Claude

Requires [uv](https://docs.astral.sh/uv/getting-started/installation/), which fetches
Python and everything else on first run.

**Claude Code:**

```bash
claude mcp add pubchem-finder -- uvx pubchem-finder mcp
```

**Claude Desktop** (Settings, Developer, Edit Config, then restart Claude):

```json
{
  "mcpServers": {
    "pubchem-finder": {
      "command": "uvx",
      "args": ["pubchem-finder", "mcp"]
    }
  }
}
```

Then ask things like:

- "What is the molecular weight and XLogP of caffeine? Is it Lipinski-compliant?"
- "Is ibuprofen hazardous? What does its GHS label say? Cite the sources."
- "What are the main targets of celecoxib in PubChem bioassays?"
- "Show me the closest analogues of aspirin and highlight what changed."
- "Does this compound have any PAINS alerts? `O=C1C(=Cc2ccccc2)SC(=S)N1`"
- With a local index: "Find the 20 compounds most similar to ibuprofen" or "Which indexed
  compounds contain an aryl sulfonamide? Draw a few with the group highlighted."

![Sulfonamide drugs with the aryl sulfonamide SMARTS match highlighted](https://raw.githubusercontent.com/MysterionRise/pubchem-finder/master/docs/images/sulfonamides-smarts.png)

*`draw_molecules` with `highlight="smarts"`: eight sulfonamide drugs with the match for
`[NX3][SX4](=O)(=O)c` highlighted, as Claude shows it inline.*

The local-index tools need an index; build one from the command line (below) and the
server picks it up. To use a non-default data directory, add
`"env": {"PUBCHEM_FINDER_HOME": "/path"}` to the Claude Desktop entry, or
`--env PUBCHEM_FINDER_HOME=/path` to `claude mcp add`.

### Tools

| Tool | What it does |
|---|---|
| `pubchem_lookup(query, by, limit)` | Resolve a name, CID, SMILES or InChIKey to PubChem records: identifiers, formula, MW, XLogP, H-bond counts, TPSA, link. |
| `compound_properties(smiles)` | RDKit descriptors and Lipinski check for any SMILES. Offline. |
| `standardize(smiles)` | Remove salts and counter-ions, neutralize, canonicalize. Offline. |
| `compare(smiles_a, smiles_b)` | Tanimoto similarity on ECFP4 and MACCS fingerprints. Offline. |
| `medchem_profile(smiles)` | QED, synthetic accessibility, Lipinski/Veber, PAINS/Brenk/NIH alerts with atoms, Murcko scaffold. Offline. |
| `draw_molecules(smiles, legends, highlight, smarts, reference_smiles)` | Structures as an inline image, highlighting a SMARTS match, what changed relative to a reference, or a similarity map. Offline. |
| `similarity_search(smiles, threshold, top_k, metric, group_by_scaffold)` | Most similar compounds in the local index, optionally grouped by Murcko scaffold. |
| `substructure_search(query, query_type, max_hits)` | Indexed compounds containing a SMILES or SMARTS fragment, lowest CID first. |
| `index_status()` | Whether the local index is built, what it holds, and how to build what is missing. |
| `pubchem_similarity_search(smiles, threshold, max_records)` | Similar compounds across all of PubChem, live, re-scored with ECFP4. |
| `pubchem_substructure_search(smiles, max_records)` | Compounds across all of PubChem that contain a fragment, live. |
| `pubchem_exact_match(smiles, identity_type)` | The PubChem record(s) for exactly this structure, optionally ignoring stereo or isotopes. |
| `pubchem_synonyms(cid, limit)` | Common, brand and systematic names, plus CAS Registry Numbers. |
| `pubchem_xrefs(cid)` | Registry IDs (CAS, ChEMBL, DrugBank...), PubMed articles and patents. |
| `pubchem_description(cid)` | Curated descriptions, with sources. |
| `pubchem_annotations(cid, heading, max_chars)` | One section of the compound page: GHS classification, hazards, pharmacology, mechanism of action, drug information, experimental properties, toxicity, use, food additives. Every item names its source. |
| `pubchem_bioactivity(cid)` | Bioassay summary: assay count, outcomes, top targets, most potent active results. |

## Use it from the command line

```bash
uvx pubchem-finder lookup caffeine                     # PubChem record by name, CID, SMILES or InChIKey
uvx pubchem-finder props "CC(=O)Oc1ccccc1C(=O)O"       # RDKit descriptors (offline)
```

`uvx` runs the latest release without installing it. To install the command once:
`uv tool install pubchem-finder` (or `pip install pubchem-finder`). The examples below
assume it is installed. To try the development version, replace `pubchem-finder` with
`--from git+https://github.com/MysterionRise/pubchem-finder pubchem-finder` in any `uvx`
command.

### Build a local index and search it

```bash
pubchem-finder index build --limit 1 --with-substructure   # first PubChem chunk: ~445,000 compounds
pubchem-finder search "CC(=O)Oc1ccccc1C(=O)O" --threshold 0.7 --top 20
pubchem-finder search "CC(=O)Oc1ccccc1C(=O)O" --group-by-scaffold
pubchem-finder substructure "c1ccc2[nH]ccc2c1"              # everything indexed that contains indole
pubchem-finder substructure "[NX3][SX4](=O)(=O)c" --smarts  # aryl sulfonamides, as SMARTS
pubchem-finder index status

pubchem-finder index build --sdf my_library.sdf.gz         # or index your own SDF files
```

`index build --limit N` downloads the first N of PubChem's ~350 SDF chunks (about 500,000
CIDs each) from NCBI over HTTPS, verifies each file's MD5, fingerprints it, deletes the
download and merges everything into one searchable index. Re-running is safe: finished
chunks are skipped, so an interrupted build resumes and a larger `--limit` extends an
existing index. `--all` indexes all of PubChem. The default fingerprint is ECFP4 (Morgan
radius 2, 1024 bits); `--radius` and `--fp-size` change it for a new index.

Substructure search is optional because it needs a second fingerprint database. It is
built from the SMILES already in the index, so `index build --with-substructure` adds it to
an existing index without downloading anything, and later builds keep it up to date
(`--without-substructure` removes it). Each compound's RDKit pattern fingerprint is
screened first, then candidates are matched atom by atom (stereo is ignored). Hits come
lowest CID first, which puts long-known compounds first; the search stops after `--max`
hits or `--max-candidates` checked compounds and says whether it saw every candidate.

### Ask PubChem directly

The `pubchem` commands query PubChem's live services, with no index needed:

```bash
pubchem-finder pubchem similarity "CC(=O)Oc1ccccc1C(=O)O" --threshold 90   # across all of PubChem
pubchem-finder pubchem substructure "c1ccc2[nH]ccc2c1" --max 50
pubchem-finder pubchem exact "CC(=O)Oc1ccccc1C(=O)O" --identity same_connectivity
pubchem-finder pubchem synonyms 2244          # names and CAS numbers
pubchem-finder pubchem xrefs 2244             # registry IDs, PubMed articles, patents
pubchem-finder pubchem describe 2244          # curated descriptions with sources
pubchem-finder pubchem annotations 2244 --heading "GHS Classification"
pubchem-finder pubchem bioactivity 2244       # assay outcomes, top targets, most potent results
```

PubChem's similarity search uses its own fingerprint and does not report scores, so each
hit is scored here with ECFP4, the same fingerprint as the local index, and sorted by it.
Successful responses are cached for 7 days, so repeated questions cost nothing and stay
within PubChem's rate limits (`pubchem-finder cache info | clear`; `PUBCHEM_FINDER_CACHE=0`
turns the cache off).

### Pictures, triage and whole files

```bash
pubchem-finder draw "CC(=O)Oc1ccccc1C(=O)O" "COC(=O)c1ccccc1O" -l aspirin -l "methyl salicylate" \
    --highlight diff --reference "CC(=O)Oc1ccccc1C(=O)O" -o analogues.png
pubchem-finder draw "O=C(O)c1ccccc1O" "CC(=O)Nc1ccc(O)cc1" --highlight similarity_map \
    --reference "CC(=O)Oc1ccccc1C(=O)O" -o similarity.png
pubchem-finder triage "O=C1C(=Cc2ccccc2)SC(=S)N1"     # QED, SA score, rules, PAINS/Brenk/NIH alerts

# .smi ("SMILES [id]" per line), .csv/.tsv with a smiles column, or - for stdin; CSV out
pubchem-finder triage -i library.smi -o triage.csv
pubchem-finder props -i library.csv -o props.csv
pubchem-finder search -i queries.smi -o hits.csv --threshold 0.6 --top 50
pubchem-finder substructure -i fragments.smi -o hits.csv
```

Rows that fail in batch mode get an `error` column instead of stopping the run.

![Aspirin analogues, with the atoms each one adds to aspirin highlighted](https://raw.githubusercontent.com/MysterionRise/pubchem-finder/master/docs/images/analogues-diff.png)

*`--highlight diff`: each analogue is aligned to the core it shares with aspirin, and the
atoms it adds are highlighted.*

![Similarity maps of four compounds against aspirin](https://raw.githubusercontent.com/MysterionRise/pubchem-finder/master/docs/images/similarity-map.png)

*`--highlight similarity_map`: atoms that raise ECFP4 similarity to aspirin in green,
atoms that lower it in pink.*

## How big is it?

The first column comes from real runs of `index build --limit 1`: on a MacBook Pro, and
with `--with-substructure` in the project's
[real-data workflow](https://github.com/MysterionRise/pubchem-finder/blob/master/.github/workflows/e2e-real.yml) on a standard GitHub Actions Linux
runner, which then runs every command and MCP tool against the result. The full-PubChem column
scales those numbers to PubChem's roughly 120 million compounds and is an estimate.

| | First PubChem chunk (measured) | All of PubChem (estimated) |
|---|---|---|
| Compounds indexed | 444,794 (29 skipped by RDKit) | ~120 million |
| Fingerprinting | 80 s on one MacBook core, 169 s on a GitHub runner | ~1 h on 8 cores, plus download time |
| Download | one SDF file, deleted after indexing | ~100 GB over HTTPS, streamed chunk by chunk |
| Similarity index file | 36 MB | ~10 GB |
| RAM for in-memory similarity search | 64 MB | ~17 GB |
| Substructure search (optional) | 169 s more to build, 88 MB file, 121 MB RAM | ~1.5 h on 8 cores, ~24 GB file, ~33 GB RAM |
| Free disk while building | small | ~45 GB, plus ~50 GB with substructure search |

Query times on the first chunk (GitHub runner): the similarity scan takes about 3 ms in
memory and 0.12 s on disk, and fetching names and SMILES for 40 hits brings a full call to
about 60 ms; a 50-hit substructure search takes 35-50 ms in memory and 0.3 s on disk. The first query also loads the index, which takes well under a second at this size.
Screening time grows roughly linearly with the size of the index.

The index is searched in memory when it fits comfortably in available RAM and from disk
otherwise; set `PUBCHEM_FINDER_INDEX_MODE=memory|disk` to force a mode. On a laptop, index
a subset with `--limit` rather than all of PubChem.

Data lives in the platform data directory (`~/.local/share/pubchem-finder` on Linux,
`~/Library/Application Support/pubchem-finder` on macOS, `%LOCALAPPDATA%\pubchem-finder`
on Windows). Override it with `--home DIR` or `PUBCHEM_FINDER_HOME`.

## Use it from Python

The CLI and MCP server are thin layers over `pubchem_finder.service`:

```python
from pubchem_finder import service

service.properties("CC(=O)Oc1ccccc1C(=O)O").molecular_weight  # 180.16
service.standardize("CC(=O)[O-].[Na+]").smiles  # 'CC(=O)O'
service.compare("CCO", "CCCO").ecfp4_tanimoto  # 0.556
service.lookup("caffeine").compounds[0].cid  # 2519 (network)
service.similarity_search("CC(=O)Oc1ccccc1C(=O)O", threshold=0.7)  # needs an index
service.substructure_search("[NX3][SX4](=O)(=O)c", query_type="smarts")
```

## Troubleshooting

- **The first run takes a minute.** uv downloads Python if needed and a few hundred MB of
  packages (RDKit, PyArrow, DuckDB...) once; later runs start in about a second.
- **Claude Desktop says the server failed to start.** Desktop apps on macOS do not see
  your shell's `PATH`, so `uvx` may not be found. Use its full path as the `command`
  (run `which uvx`, typically `~/.local/bin/uvx`). The server's log is in Claude
  Desktop's MCP logs.
- **"No local index"** from `similarity_search` or `substructure_search`: build one with
  `pubchem-finder index build --limit 1` (add `--with-substructure` for substructure
  search). `index_status` says what is missing.
- **PubChem is slow or busy.** Requests are throttled to PubChem's limit of 5 per second
  and retried when PubChem is busy. Structure searches across all of PubChem can time out
  on very generic fragments; use a more specific query, or the local index.
- **Behind a proxy:** set `HTTPS_PROXY`; all requests to PubChem and NCBI honour it.
- **Linux on ARM (aarch64):** FPSim2 publishes no wheel there, so installing needs a C++
  compiler to build it.
- **Starting over:** `pubchem-finder cache clear` empties the response cache; deleting the
  data directory removes the index too.

## Data, credits and citation

The data comes from PubChem, the National Library of Medicine's open chemistry database.
Please follow [NCBI's policies](https://www.ncbi.nlm.nih.gov/home/about/policies/) and
[PubChem's guidance on citing](https://pubchem.ncbi.nlm.nih.gov/docs/citation-guidelines),
and cite PubChem when you use its data:

> Kim S, et al. PubChem 2025 update. *Nucleic Acids Research* 2025;53(D1):D1516-D1525.
> [doi:10.1093/nar/gkae1059](https://doi.org/10.1093/nar/gkae1059)

Annotations (safety, pharmacology, toxicity...) are deposited by other organisations;
each one names its source, and some carry their own license terms.

PubChem Finder is built on [RDKit](https://www.rdkit.org/),
[FPSim2](https://github.com/chembl/FPSim2), [DuckDB](https://duckdb.org/),
[Apache Arrow](https://arrow.apache.org/) and the
[MCP Python SDK](https://github.com/modelcontextprotocol/python-sdk). To cite PubChem
Finder itself, see [CITATION.cff](https://github.com/MysterionRise/pubchem-finder/blob/master/CITATION.cff) (GitHub's "Cite this repository" button).

## Development

Requires [uv](https://docs.astral.sh/uv/). Python 3.11, 3.12 and 3.13 are supported; the
repository's `.python-version` makes `uv sync` use 3.12 even if pyenv or conda put an
older Python first on your PATH.

```bash
uv sync            # create .venv and install everything, including dev tools
make check         # ruff lint + format check, mypy (strict), pytest with coverage
make test          # just the tests
PUBCHEM_LIVE=1 uv run pytest -m live   # tests that call the real PubChem services
uv run python scripts/smoke.py --live  # every MCP tool through the real stdio server
```

The animation at the top is drawn from a real index by `scripts/demo_media.py`
(`collect` runs the searches and times them, `render` writes the GIF and MP4); its
docstring has the exact commands, and the real-data workflow renders it on every run.
`docs/media/demo-data.json` holds the results behind the current animation, so
`uv run python scripts/demo_media.py render docs/media/demo-data.json -o docs/media`
redraws it without an index.

CI runs the tests on Linux, macOS and Windows, smoke-tests the built wheel on all three,
and runs the live PubChem tests weekly. See [CONTRIBUTING](https://github.com/MysterionRise/pubchem-finder/blob/master/.github/CONTRIBUTING.md),
[ROADMAP](https://github.com/MysterionRise/pubchem-finder/blob/master/ROADMAP.md) and [CHANGELOG](https://github.com/MysterionRise/pubchem-finder/blob/master/CHANGELOG.md).

```
PubChem Compound SDF (HTTPS) --> streaming RDKit parse --> chunks/<chunk>.parquet  (CID, SMILES, InChIKey, name, formula, MW)
                                                      \-> chunks/<chunk>.h5       (FPSim2 fingerprints)
                                                            \-> merge + popcount sort -> index/pubchem.h5 + manifest.json
chunks/<chunk>.parquet SMILES --> RDKit pattern fingerprints --> chunks/pattern/<chunk>.h5 -> index/pubchem-pattern.h5
                                                                  (only with --with-substructure)
pubchem-finder CLI  ----\                 /-> FPSim2 exact similarity (in memory or on disk) + DuckDB over Parquet
                         +--> service ---+--> FPSim2 pattern screen + RDKit atom match (substructure)
MCP server (stdio) -----/                 \-> PubChem PUG REST / PUG View for live data
```

This repository started in 2020 as an Indigo + Elasticsearch indexer for PubChem dumps;
that code is in the git history before commit `06619d9`.

## License

MIT. See [LICENSE](https://github.com/MysterionRise/pubchem-finder/blob/master/LICENSE).
