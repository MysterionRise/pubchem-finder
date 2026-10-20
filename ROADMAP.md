# Roadmap

Decisions: FPSim2 + DuckDB/Parquet as the search backend (no server required),
MCP over stdio, uv + hatchling + ruff, Python 3.11 to 3.13.

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

## Phase 4: PubChem extensions (live PUG REST / PUG View, no local data)

- [ ] Remote structure search: `pubchem_similarity_search`, `pubchem_substructure_search`,
      `pubchem_exact_match` (PubChem's fast search over all of PubChem)
- [ ] `pubchem_synonyms` (names, CAS numbers) and `pubchem_xrefs` (patents, PubMed, registry IDs)
- [ ] `pubchem_description` and `pubchem_annotations` for curated PUG View sections
      (GHS classification, pharmacology, drug information, experimental properties, toxicity)
- [ ] `pubchem_bioactivity` (aggregated assay summary) and `pubchem_structure_image` (PNG)
- [ ] Matching `pubchem-finder pubchem ...` CLI commands

## Phase 5 (optional): PostgreSQL + RDKit cartridge backend

- [ ] `SearchBackend` implementation with substructure search, docker-compose,
      `substructure_search` MCP tool
