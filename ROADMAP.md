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

## Phase 2: PUG REST + MCP server (v0.1, works with no local index)

- [ ] `data/pug_rest.py`: PubChem PUG REST client (name / CID / SMILES / InChIKey)
- [ ] `mcp_server.py` (stdio) with `pubchem_lookup`, `compound_properties`,
      `standardize`, `index_status`, `similarity_search` (errors until an index exists)
- [ ] `pubchem-finder lookup | props | standardize | mcp` CLI
- [ ] README: `claude mcp add` and Claude Desktop configuration

## Phase 3: local index + similarity search (v0.2)

- [ ] HTTPS download of `Compound/CURRENT-Full/SDF/*.sdf.gz` with resume and md5 check
- [ ] Streaming SDF parse, per-chunk Parquet metadata + FPSim2 `.h5`, `merge_db_files`
- [ ] `index build | finalize | status`, `search --smiles ... --threshold ... --top ...`
- [ ] In-memory vs on-disk FPSim2 engine selection; lazy load in the MCP server
- [ ] Honest sizing in the README (full PubChem is an overnight job, ~20 GB RAM in-memory)

## Phase 4 (optional): PostgreSQL + RDKit cartridge backend

- [ ] `SearchBackend` implementation with substructure search, docker-compose,
      `substructure_search` MCP tool
