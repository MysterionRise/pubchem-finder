# Roadmap

Decisions: FPSim2 + DuckDB/Parquet as the search backend (no server required),
MCP over stdio (streamable HTTP optional), uv + hatchling + ruff, Python 3.11 to 3.14.

## Phase 0: repo hygiene (done)

- [x] Remove the 2020 Indigo/Elasticsearch code, cookiecutter template and stale docs
- [x] uv + hatchling project, ruff, strict mypy, single CI workflow on 3.11/3.12/3.13
- [x] Recover the integration tests that were never collected

## Phase 1: chemistry core (done)

- [x] Morgan fingerprints via `rdFingerprintGenerator` (parity test against RDKit)
- [x] Standardizer order: cleanup, largest organic fragment, then neutralize
- [x] Vectorized Tanimoto/Dice, cosine, RDKit bulk Tanimoto helper
- [x] `Molecule` with identity equality and timezone-aware timestamps

## Phase 2: PUG REST + MCP server (v0.1, works with no local index) (done)

- [x] `data/pug_rest.py`: PubChem PUG REST client (name / CID / SMILES / InChIKey),
      throttled to 5 req/s, retries on 429/503, typed errors
- [x] `mcp_server.py` (stdio) with `pubchem_lookup`, `compound_properties`,
      `standardize`, `compare`
- [x] `pubchem-finder lookup | props | standardize | compare | mcp` CLI
- [x] README: `uvx` quick start, `claude mcp add` and Claude Desktop configuration

## Phase 3: local index + similarity search (v0.2) (done)

- [x] HTTPS download of `Compound/CURRENT-Full/SDF/*.sdf.gz` with resume and MD5 check
- [x] Streaming SDF parse; per-chunk Parquet metadata + FPSim2 `.h5`; parallel, resumable build
- [x] Merge and popcount-sort into one index, skipped when already up to date
- [x] `index build | status` and `search` CLI; `--sdf` to index your own files
- [x] In-memory vs on-disk engine selection; lazy load in the MCP server
- [x] MCP tools `similarity_search` and `index_status`
- [x] Measured throughput and sizing in the README

## Phase 4a: PubChem knowledge tools (done)

- [x] Remote structure search: `pubchem_similarity_search` (re-scored with ECFP4),
      `pubchem_substructure_search`, `pubchem_exact_match`
- [x] `pubchem_synonyms` (names, validated CAS numbers) and `pubchem_xrefs`
      (registry IDs, PubMed, patents)
- [x] `pubchem_description` and `pubchem_annotations` for curated PUG View sections
      (GHS classification, pharmacology, drug information, experimental properties,
      toxicity...), each annotation attributed to its source
- [x] `pubchem_bioactivity`: assay outcomes, top targets, most potent results
- [x] 7-day on-disk response cache; `cache info | clear`
- [x] `pubchem-finder pubchem ...` CLI commands

## Phase 4b: assistant features (done)

- [x] `draw_molecules`: RDKit structure grids as inline images, with SMARTS, "what changed"
      (MCS, aligned to the shared core) and similarity-map highlighting; `pubchem-finder draw`
- [x] `medchem_profile`: QED, synthetic accessibility, Lipinski/Veber, PAINS/Brenk/NIH
      alerts, Murcko scaffold; scaffold grouping of similarity hits; `pubchem-finder triage`
- [x] Batch CSV mode for `props`, `standardize`, `triage` and `search`

## Phase 5: local substructure search (done)

- [x] Optional RDKit pattern-fingerprint database (`index build --with-substructure`), built
      from the stored SMILES, so it can be added to an existing index offline
- [x] FPSim2 screen, then RDKit verification in CID order, with hit and candidate caps
- [x] `substructure_search` MCP tool and `pubchem-finder substructure` CLI (SMILES or SMARTS,
      batch files)
- [x] Faster index finalize: popcount bins computed with numpy instead of FPSim2's per-row loop

## Phase 6: tested everywhere, released as 2.0.0 (done)

- [x] Test suite on Linux, macOS and Windows; the built wheel installed and smoke-tested
      on all three; the real MCP server tested over stdio
- [x] Live PubChem tests weekly; an on-demand workflow that builds an index from real
      PubChem files and calls every CLI command and MCP tool
- [x] Memory detection on macOS and Windows (psutil), readable atom labels in drawings
- [x] PyPI release through trusted publishing, CHANGELOG, citation file, issue templates,
      security policy, Dependabot

## Phase 7: listed, hardened, 2.1 (done)

- [x] Streamable HTTP transport (#33): `pubchem-finder mcp --transport http`
- [x] Listed in the MCP registry (`server.json`, published by the release workflow) and a
      Claude plugin (`plugins/pubchem-finder`) for Claude Code and Cowork, submitted to
      Anthropic's plugin directory
- [x] 2.1.1: tool titles and annotations (read-only, open-world), `--python 3.14` in every
      published launch command (Python 3.15 lacks FPSim2, numexpr and ndindex wheels), an
      index that survives RDKit/FPSim2 upgrades, each CID once from overlapping chunks,
      readable errors and a cache that falls back to the network
- [x] Release workflow guarded (master only, CI-green commits, tag before PyPI, approval
      for PyPI), `make release-prep`, the plugin bump prepared by each release
      ([docs/releasing.md](docs/releasing.md))
- [x] CI gates: SHA-pinned actions, zizmor and actionlint, ruff security rules, a 95%
      coverage floor (branches and server subprocesses), the lowest allowed dependency
      versions, a monthly real-data run, issues for failing scheduled runs, CodeQL
- [x] Tested VS Code (GitHub Copilot) setup in [docs/clients.md](docs/clients.md), the
      first from an outside contributor (#37)

## 2.2: the first five minutes (done)

- [x] A demo index of 1,932 FDA-approved drugs shipped with the package: local
      similarity and substructure search work on first use, in about a second, and say
      they searched the demo set until a PubChem index is built (`index build --demo`)
- [x] Hints that run where users are: `uvx --python 3.14 pubchem-finder==<version> index
      build ...` with honest sizes (PubChem's first chunk is a 323 MB download), and
      download progress
- [x] MCP prompts (compound dossier, find analogues, triage) and resources
      (`pubchem://compound/{cid}`, index status)
- [x] A format version in the index manifest, checked on read

## 2.2.1 and 2.2.2: right answers (done)

- [x] 2.2.1: in-memory searches of one index run one at a time (a race in FPSim2,
      chembl/FPSim2#140); size limits on structures, so an oversized SMILES is an error
      and not a crash; batch mode cannot overwrite its input; `pubchem_xrefs` fetches
      patents on their own
- [x] 2.2.2: `--sdf` rebuilds that survive failures (the index records which build of
      each chunk it holds, and nothing is deleted before its replacement is complete);
      substructure search kept through an interrupted merge; bioassay summaries that rank
      only potencies and count each gene once; a CID PubChem does not have is "not
      found"; record IDs the index cannot hold are skipped

## 2.3: first run (next)

- [ ] One-click install for Claude Desktop: an `.mcpb` bundle attached to each release
      (#109)
- [ ] Ready-made indexes: `index fetch` downloads a prebuilt slice of PubChem instead of
      building one from up to 111 GB (#110)
- [ ] Damaged index files named in the error and repaired by `index build` (#111); first
      searches from two programs at once (#112)
- [ ] Batch mode reads CSV files saved from Excel (#113); one-line errors for unwritable
      paths (#114)

## 2.4: from a list to an answer

- [ ] Resolve a list of names, CAS numbers or InChIKeys in one call, with what did not
      resolve (#115)
- [ ] Check structures an assistant proposed: valid, canonical, no duplicates (#116)
- [ ] SAR tools: clusters, a diverse subset, matched molecular pairs, an R-group table
      (#117)

## 2.5: see the results

- [ ] Search hits as an interactive grid in the chat, in clients that support MCP Apps
      (#118)

## 2.6: predictions

- [ ] ADMET predictions as an optional extra, each value marked as a prediction (#119)

## Ongoing

- [ ] Tested setups for more MCP clients (Cursor, Windsurf, Zed, JetBrains, Antigravity /
      Gemini CLI, Codex CLI, Cline, LM Studio): #23, one good first issue per client
- [ ] Examples for agent frameworks: LangChain (#34), OpenAI Agents SDK (#106),
      PydanticAI or LlamaIndex (#107)

## Kept, not scheduled

- Upkeep of indexes you build yourself: retry failed downloads (#121), expiry and a size
  cap for the response cache (#122), memory of very broad substructure queries (#123)
- Python 3.15, once FPSim2 publishes wheels (#124; a weekly job watches)
- Shared servers: a Docker image, `--allowed-host`, per-request limits and time budgets.
  Designed, and scheduled when a team asks: #125

## Ideas with evidence behind them

- Bioactivity by target ("what is active against EGFR"), probably from ChEMBL
- Agent Skills for the common workflows (dossier, analogue search, triage)
- 3D structures: PubChem's conformers, or RDKit-generated ones
- A try-it-in-the-browser demo of the built-in drug set
- A short paper (JOSS or ChemRxiv)
