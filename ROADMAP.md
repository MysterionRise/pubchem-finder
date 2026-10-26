# Roadmap

Decisions: FPSim2 + DuckDB/Parquet as the search backend (no server required),
MCP over stdio (streamable HTTP optional), uv + hatchling + ruff, Python 3.11 to 3.13.

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

## Next, driven by feedback

- [ ] Tested setups for MCP clients beyond Claude (Cursor, VS Code, Windsurf, Zed,
      JetBrains, Antigravity / Gemini CLI, Codex CLI, Cline, LM Studio) in
      [docs/clients.md](docs/clients.md): #23, one good first issue per client
- [x] List the server in the MCP registry (`server.json`, published by the release
      workflow) and ship a Claude plugin (`plugins/pubchem-finder`) for Claude Code and Cowork
- [ ] A `.mcpb` bundle for one-click install in Claude Desktop, if people ask: Anthropic's
      directory no longer lists `.mcpb` files, so it would be a release asset only
- [ ] A quick demo index (a few thousand well-known compounds via PUG REST), so the first
      local search takes seconds instead of a 250 MB download
- [ ] MCP prompts and resources: a compound dossier prompt, `pubchem://compound/{cid}`
- [ ] Matched molecular pairs across a hit set (`rdMMPA`)
- [ ] ADMET predictions as an optional extra (pulls in PyTorch)
- [x] Streamable HTTP transport (#33): `pubchem-finder mcp --transport http`
- [ ] A Docker image, and `--allowed-host` for the HTTP transport, if a shared server is wanted
