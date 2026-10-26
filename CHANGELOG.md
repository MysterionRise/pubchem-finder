# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project uses
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Changed

- The Claude plugin runs the released package from its own lockfile
  (`plugins/pubchem-finder/pyproject.toml` and `uv.lock`: every dependency at an exact
  version with a checked hash) and has an icon and documentation and support links, as
  Anthropic's plugin directory requires. `make plugin-lock VERSION=...` points it at a new
  release; tests keep its files in step and never ahead of the package.

## [2.1.0] - 2026-10-26

### Added

- **MCP over streamable HTTP:** `pubchem-finder mcp --transport http [--host 127.0.0.1]
  [--port 8000]` serves the same 17 tools at `http://HOST:PORT/mcp` for clients that
  connect to a URL and for one shared index on a team machine. stdio stays the default.
  `GET /health` reports the version. The server has no authentication or TLS of its own;
  docs/clients.md, "Over HTTP", explains the reverse-proxy setup and the shared PubChem
  rate limit ([#33](https://github.com/MysterionRise/pubchem-finder/issues/33)).
- `pubchem-finder mcp --home DIR`, like the other commands.
- The server is listed in the [MCP registry](https://registry.modelcontextprotocol.io)
  as `io.github.MysterionRise/pubchem-finder` (`server.json`; the release workflow
  publishes each version), and a Claude plugin (`/plugin marketplace add
  MysterionRise/pubchem-finder`) installs it in Claude Code and Cowork.
- Python 3.14 is tested in CI and listed in the classifiers.
- `examples/mcp_client.py`: the MCP tools from plain Python with the official client, the
  starting point for agent-framework examples ([#34](https://github.com/MysterionRise/pubchem-finder/issues/34)).
- Reviewing client-setup contributions: the pull request template asks for the tool
  calls as the client showed them, docs/clients.md says what reviewers check, and a
  repository skill (`.claude/skills/review-client-setup`) reproduces the reported
  answers with the project before the maintainer posts a review.
- `scripts/smoke.py --url` runs the MCP checks against a running HTTP server; CI does so
  with the installed wheel, and an end-to-end test starts the HTTP server on Linux, macOS
  and Windows.

### Changed

- PubChem's similarity and substructure hits are now re-scored with the local index's
  default fingerprint (ECFP4, 1024 bits) instead of 2048 bits, so `ecfp4_similarity` is
  comparable with `similarity_search` as the docs said. `compare` keeps 2048 bits and says so.

- The unused `slow` test marker is gone, so `pytest` runs everything (the suite takes
  about ten seconds). The source distribution no longer ships the demo media.

### Fixed

- CITATION.cff said 2.0.0; a test now keeps it, `server.json` and the plugin manifest
  in step with the package version.
- Concurrent MCP tool calls could exceed PubChem's limit of 5 requests per second: the
  shared client's throttle now holds a lock across its wait.

## [2.0.1] - 2026-10-22

Documentation release: the package's code and behaviour are unchanged.

### Added

- An animated walkthrough at the top of the README: a similarity search for aspirin and a
  SMARTS substructure search, rendered from real results on PubChem's first 444,794
  compounds. `scripts/demo_media.py` collects the results from an index and renders the
  GIF and MP4 (`docs/media/`); the real-data workflow renders it on every run.
- `docs/clients.md`: how to run the MCP server in clients other than Claude (Cursor,
  VS Code, Windsurf, Zed, JetBrains, Codex, LM Studio...), with a quick check and a
  status table. Tested setups are welcome, see
  [#23](https://github.com/MysterionRise/pubchem-finder/issues/23).

### Fixed

- Images and links in the README now use absolute URLs, so they show on PyPI too.
- A download test that failed in 1 run in 256.

## [2.0.0] - 2026-10-22

A complete rewrite. PubChem Finder is now a command-line tool, a Python library and an
MCP server that lets Claude (or any MCP client) search PubChem by chemical structure,
with a local index you can build on a laptop.

### Added

- **MCP server over stdio** (`pubchem-finder mcp`) with 17 tools for Claude Desktop,
  Claude Code and other MCP clients.
- **Local similarity search.** `index build` downloads PubChem's SDF files (or indexes
  your own), fingerprints them with RDKit (ECFP4 by default) into an
  [FPSim2](https://github.com/chembl/FPSim2) index with Parquet metadata, and searches
  it exactly, in memory or from disk. Builds run in parallel and resume after an
  interruption.
- **Local substructure search** (`index build --with-substructure`, `substructure`):
  an RDKit pattern-fingerprint screen, then an exact atom-by-atom match. SMILES or
  SMARTS queries, results in CID order.
- **PubChem knowledge tools:** lookup by name, CID, SMILES or InChIKey; PubChem's own
  similarity, substructure and identity searches (hits re-scored with ECFP4); synonyms
  and CAS numbers; cross-references (registry IDs, PubMed, patents); descriptions;
  curated page sections such as GHS classification, pharmacology and toxicity, each
  attributed to its source; bioassay summaries. Responses are cached for 7 days.
- **Chemistry tools:** RDKit descriptors and Lipinski check, standardization, ECFP4 and
  MACCS comparison, med-chem triage (QED, synthetic accessibility, Lipinski/Veber,
  PAINS/Brenk/NIH alerts, Murcko scaffold), and structure drawings that highlight a
  SMARTS pattern, what changed relative to a reference, or a similarity map.
- **Batch mode** for `props`, `standardize`, `triage`, `search` and `substructure`:
  `.smi`, `.csv` or `.tsv` in, CSV out, with an error column for bad rows.
- **Testing:** the suite runs on Linux, macOS and Windows; the built wheel is installed
  and smoke-tested on all three; live PubChem tests run weekly; a real-data workflow
  builds an index from PubChem's own files and calls every tool.

### Changed

- Python 3.11 to 3.13, packaged with uv and hatchling, linted with ruff, type-checked
  with mypy in strict mode.

### Removed

- The 2020 Indigo + Elasticsearch indexer, the cookiecutter template and the Docker
  files. They remain in the git history before commit `06619d9`.

## 1.x (2020)

An Indigo + Elasticsearch indexer for PubChem SDF dumps. It was never released as a
package.

[Unreleased]: https://github.com/MysterionRise/pubchem-finder/compare/v2.1.0...HEAD
[2.1.0]: https://github.com/MysterionRise/pubchem-finder/releases/tag/v2.1.0
[2.0.1]: https://github.com/MysterionRise/pubchem-finder/releases/tag/v2.0.1
[2.0.0]: https://github.com/MysterionRise/pubchem-finder/releases/tag/v2.0.0
