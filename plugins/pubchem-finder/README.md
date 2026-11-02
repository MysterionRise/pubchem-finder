# PubChem Finder

Search [PubChem](https://pubchem.ncbi.nlm.nih.gov/) by chemical structure from Claude.
This plugin starts the [pubchem-finder](https://github.com/MysterionRise/pubchem-finder)
MCP server, which gives Claude 17 tools: PubChem lookups by name, CID, SMILES or InChIKey;
synonyms and CAS numbers; GHS safety, pharmacology and toxicity sections with their
sources; bioassay summaries; PubChem's own similarity, substructure and exact-match
searches; RDKit descriptors, standardization, med-chem triage and structure drawings;
and exact local similarity and substructure search with no rate limits, over a built-in
set of FDA-approved drugs from the start and over PubChem's compounds once you build an
index.

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
one pinned release of the package (`uvx --python=3.14 pubchem-finder==<version> mcp`);
the first start downloads it and its dependencies (RDKit, PyArrow, DuckDB) and takes
about a minute, later starts a few seconds. `--python=3.14` keeps uv off Python 3.15, for
which FPSim2, numexpr and ndindex publish no wheels yet. The lockfile next to the plugin
records the exact dependency set Anthropic's directory verified. Then ask things like:

- "What is the molecular weight and XLogP of caffeine? Is it Lipinski-compliant?"
- "Is ibuprofen hazardous? What does its GHS label say? Cite the sources."
- "Show me the closest analogues of aspirin and highlight what changed."

Local similarity and substructure search work from the first question: until you build
an index they search a built-in demo set of about 1,900 FDA-approved drugs (set up in a
few seconds on first use), and each answer says so. To search PubChem's compounds
locally, build an index once from a terminal:
`uvx --python 3.14 pubchem-finder index build --limit 1 --with-substructure` (PubChem's
first 445,000 compounds: a 323 MB download and about 5 minutes; `--all` for
everything). The server switches to it automatically.

The server also offers three prompts in the `/` menu (`compound_dossier`,
`find_analogues`, `triage_structures`) and PubChem compound records as resources
(`pubchem://compound/<cid>`) that you can attach with `@`.

## Data

The server runs on your machine. The `pubchem_*` tools send your queries (names,
identifiers, structures) to PubChem's public web services at ncbi.nlm.nih.gov and cache
the responses locally for 7 days; the RDKit tools, the demo set and the local index work
offline.
Nothing is sent anywhere else.
