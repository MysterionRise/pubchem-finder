# PubChem Finder

Search [PubChem](https://pubchem.ncbi.nlm.nih.gov/) by chemical structure from Claude.
This plugin starts the [pubchem-finder](https://github.com/MysterionRise/pubchem-finder)
MCP server, which gives Claude 17 tools: PubChem lookups by name, CID, SMILES or InChIKey;
synonyms and CAS numbers; GHS safety, pharmacology and toxicity sections with their
sources; bioassay summaries; PubChem's own similarity, substructure and exact-match
searches; RDKit descriptors, standardization, med-chem triage and structure drawings;
and, once you build a local index, exact similarity and substructure search over it with
no rate limits.

## Install

In Claude Code (2.1.275 or later; older versions take the two-step form below):

```
/plugin install pubchem-finder --marketplace MysterionRise/pubchem-finder
```

```
/plugin marketplace add MysterionRise/pubchem-finder
/plugin install pubchem-finder@mysterionrise
```

## Use it

Requires [uv](https://docs.astral.sh/uv/getting-started/installation/). The plugin runs
the released package from its own lockfile, so the first start installs pinned,
hash-checked dependencies (RDKit, PyArrow, DuckDB) and takes about a minute; later starts
take under a second. It uses Python 3.12 to 3.14 and uv fetches one if none is installed.
Then ask things like:

- "What is the molecular weight and XLogP of caffeine? Is it Lipinski-compliant?"
- "Is ibuprofen hazardous? What does its GHS label say? Cite the sources."
- "Show me the closest analogues of aspirin and highlight what changed."

For local structure search, build an index once from a terminal:
`uvx pubchem-finder index build --limit 1 --with-substructure` (PubChem's first 445,000
compounds; `--all` for everything). The server picks it up automatically.

## Data

The server runs on your machine. The `pubchem_*` tools send your queries (names,
identifiers, structures) to PubChem's public web services at ncbi.nlm.nih.gov and cache
the responses locally for 7 days; the RDKit tools and the local index work offline.
Nothing is sent anywhere else.
