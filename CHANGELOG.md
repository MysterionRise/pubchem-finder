# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project uses
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Fixed

- CLI commands and MCP tools can start on Linux without RDKit’s drawing libraries. Drawing reports which libraries are needed when they are missing (#175).

## [2.2.3] - 2026-11-11

### Fixed

- On an Apple-silicon Mac, everything that needs a Bemis-Murcko scaffold crashed with
  `RuntimeError: bad any cast`: `triage`, the `medchem_profile` tool, and grouping
  search hits by scaffold. RDKit 2026.9.2, published on 2026-11-11, raises it from
  `MurckoDecompose` for every molecule in its macOS wheels (its Linux wheels are fine).
  A new install picks the newest RDKit, and `uvx` can pick it up on a later start, so
  setups that had worked were hit as well. That RDKit release is now excluded until a
  fixed one is published. To repair an install that has it, run
  `pip install -U pubchem-finder`, or start it once with `uvx --refresh`. A command
  that pins an older pubchem-finder (`pubchem-finder==2.2.2`) keeps the fault until
  the pin moves.

  One side effect, on Linux without the X11 libraries (a slim container): RDKit
  2026.9.2 was the first release that starts there. The earlier releases, which
  installs get again, need `libxrender1`, `libxext6` and `libexpat1` as before; the
  README's Troubleshooting section now says so.

## [2.2.2] - 2026-11-11

### Fixed

- `pubchem_bioactivity` listed measurements that are not potencies among a compound's
  "most potent" results (#103). PubChem puts every value reported in micromolar into
  one column, so enzyme rates (Vmax, Km), plasma concentrations and solubilities were
  ranked with IC50 and Ki values: for quercetin, two Vmax values came 4th and 5th. Only
  IC50, EC50, AC50, Potency, Ki, Kd, ID50, ED50, MIC and GI50 are ranked now, each
  result carries its gene ID, and the descriptions say "lowest potency-type values"
  and to check the assay. Which measures belong on that list is open in #101.
- `top_targets` counts each gene once (#103). PubChem lists a gene under a UniProt
  accession, NCBI protein IDs or none, and each spelling was a separate target, so
  counts were too low and real top targets dropped out: for quercetin, ESR1 now has
  14 active assays (was 11) and CYP1A2 is back in the top five. A row that gives the
  accession without a gene is counted under the gene the same table gives for that
  accession (ChEMBL's rows often lack it): imatinib listed ABL1 twice, with 141 and
  49 assays, and now once with 190. A protein that PubChem lists under two gene IDs
  (isoniazid's InhA, calmodulin) is one target as well, and `most_potent` lists a
  result once however many rows repeat it: per gene, or per protein of a receptor
  (nicotine's top five were one binding assay, once per nAChR subunit). The accession
  shown for a gene is UniProt's whenever the table has one, also when only non-active
  rows carry it.
- A CID that PubChem does not have is "not found" (#104). `pubchem_lookup` by CID and
  `pubchem://compound/{cid}` returned an empty record, `pubchem_bioactivity` said "not
  tested" and `pubchem_xrefs` reported zero references, because PubChem answers a
  missing compound and a compound without data alike. The same holds for a lookup by
  SMILES or InChIKey of a structure PubChem does not have, which returned a record
  with CID 0 and a link to `/compound/0`: it now raises `CompoundNotFound`, so code
  that tested for `cid == 0` or an empty title has to catch that instead. A lookup of
  several CIDs at once returns the ones that exist and leaves the others out.
- `lookup` says how to look up something that is not a name. `lookup 2244` and
  `lookup <SMILES>` answered "No CID found ... the given name"; the error now adds
  `--by cid`, `--by smiles` or `--by inchikey`, and the README shows it.
- `draw_molecules` with `highlight="diff"` (and `draw --highlight diff`) no longer
  marks almost the whole molecule as changed. The common-substructure search has a
  2-second limit, and when it ran out the few atoms found so far were taken for the
  answer: dolutegravir drawn against cabotegravir had 28 of its 30 atoms highlighted,
  where 5 differ. About one analogue pair in eleven hit the limit. A search that was
  cut short is no longer used; a second method (RDKit's RASCAL, which answers for
  close analogues in milliseconds) is asked instead. It is not asked about
  macrocycles, not believed when its answer is implausibly small, and not asked when
  the first search was stopped by a busy machine without having worked (that search
  is repeated once); then, and with RDKit 2024.09, the panel has no highlights and
  its caption says "difference not computed". Stereoisomers and other molecules with
  the same atoms and bonds need no search at all. In Python,
  `drawing.diff_highlights` now returns a third list, the positions of the panels
  that were not computed.
- `highlight="diff"` no longer shows a changed ring as unchanged, or a molecule that
  shares nothing with the reference as identical to it. An aromatic ring matched a
  saturated one (a phenyl that became cyclohexyl was drawn without any highlight),
  and a molecule with no common substructure was drawn blank, like an exact copy.
  Atoms whose ring went from aromatic to saturated, or back, are now highlighted, and
  a molecule that shares nothing is highlighted whole.
- `standardize` keeps sulfoxides as `S=O`. RDKit's cleanup redraws them `[S+]-[O-]`,
  so the "standardized" omeprazole no longer matched omeprazole: 0.63 against itself
  in a similarity search, no substructure hit, a different QED and a Brenk alert for
  a charged sulfur that the drug does not have. 15 of the 1,932 demo drugs were
  affected, and the `triage_structures` prompt standardizes before it profiles.
- `standardize` no longer returns the counter-ion or a ligand as the parent. For
  pyrantel pamoate it returned pamoic acid, for diatrizoate meglumine the meglumine
  and for fosfomycin tromethamine the tromethamine (RDKit's "largest fragment" counts
  hydrogens); cisplatin came back as ammonia, carboplatin and oxaliplatin as their
  ligands, sodium chloride as HCl. The parent is now the organic component with the
  most heavy atoms, passing over about 70 common salt formers and solvents (pamoate,
  meglumine, tromethamine, tartrate, citrate, benzathine, carbonate, ethanol) while
  the input has another. A compound of a metal other than an alkali or alkaline-earth metal,
  drawn with bonds or as ions, and an input with nothing organic in it are returned
  as given with a `note`. The new `removed` field lists every fragment left out, and
  the note says when one of them is an organic component that is not a common
  counter-ion (the other half of a two-drug mixture) or a metal such as lithium. The
  CLI prints both on stderr, and batch output gains `removed` and `note` columns.
  Salts of zinc, iron, silver or bismuth are no longer stripped to the acid: they
  come back as given, with the note. In Python,
  `MoleculeStandardizer.standardize_detailed` returns the same three parts, and the
  `triage_structures` prompt tells the assistant to pass `removed` and `note` on.
- `index build --sdf` could leave an index that searched an old version of a file while
  reporting it ready (#102). 2.2.1 rebuilt the chunk of a changed file but decided whether to
  merge from what the same run had done, so a rebuild that failed part-way (a full disk,
  Ctrl+C, a damaged file) was followed by a run that merged nothing. Each chunk build
  now has an id, the index records the ids it holds, and `index build` merges whenever
  they differ. An index that 2.2.1 merged from `--sdf` files is merged once more the
  next time it is built on, which repairs one already in that state; the build says
  that this is why it merges.
- `index build --sdf` refuses files in which one ID stands for different compounds.
  A hit's structure is looked up by its ID, so two vendor files that each number their
  records from 1, a file that uses an ID twice, or a library numbered from 1 next to
  PubChem's chunks were merged without a word and then searched wrongly: aspirin came
  back with paracetamol's structure at similarity 1.000, and a benzene substructure
  search found 2 of the 4 compounds containing it. The build now stops before merging,
  names the ID, both structures and both files, and says which file to renumber. The
  same compound in two files is still merged, also when they draw it differently
  (Kekulé or aromatic, acid or anion, with or without stereo); a salt and its parent,
  or two tautomers, are different compounds here. An index built earlier from such
  files stays as it is until its next build.
- A rebuilt `--sdf` chunk is no longer read together with the fingerprints of the file
  it replaced. Its records went live when its files were swapped in, while the index
  kept the old fingerprints until the merge: after a build that failed or was stopped
  in between (and for a running server, during every rebuild), a hit carried whatever
  structure the new file had under the same ID, and `index status` said ready. With
  IDs that were renumbered between exports, a search for aspirin returned caffeine's
  structure at similarity 1.000. Searches now leave such a chunk's structures out and
  carry a `warning` saying which file to build again (batch mode prints it once on
  stderr); `index status` lists it in `unmerged_chunks`; a build that fails after
  replacing chunks merges them before it reports the failure; and a failed merge no
  longer says "The current index is unchanged" when a chunk was replaced. Hits are
  judged against the merge that produced them, so a build that finishes between a
  search and the lookup of its hits cannot pair them either, and a substructure search
  that had to leave candidates out no longer calls itself exhaustive.
- An `--sdf` file whose name contains `[`, `*` or `?` is read back by that name. The
  chunk files of `lib[1].sdf` were looked up as a pattern, which matches `lib1` and
  not `lib[1]`: its hits came back without a structure, or with the structure
  `lib1.sdf` has under the same ID.
- A chunk being rebuilt is replaced only when its new files are complete, and every
  check runs first (#102). 2.2.1 deleted the old chunk up front, so a build that was then
  refused (for example because RDKit had changed since the other chunks were built) or
  that failed had already lost it. A failed build also leaves no `.tmp` files behind.
- A file with no usable record does not replace a chunk that has some. An export that
  failed (0 bytes, a saved error page) or whose IDs were no longer numbers was
  fingerprinted into an empty chunk that took the place of the last good one: with one
  library the index stayed "ready" on fingerprints without records and every later
  build failed with "No built chunks"; with several, the library left the index
  without a word. The build now stops, keeps the chunk and says how to take the
  library out of the index if that is what was meant.
- `index build --sdf` keeps the chunk of a file that only moved (same size and
  modification time) instead of refusing it, compares chunk names without case
  (`Library.sdf` and `library.sdf` are the same files on macOS and Windows), and says
  so when one file is given twice. Chunks that already exist keep their spelling, so
  an index holding both `Library` and `library` (possible on Linux until now) still
  builds. A chunk built by 2.2.0 or earlier does not record its file and is rebuilt
  once from the file passed for it, as in 2.2.1. When the file a chunk was built from
  is gone and another one has its name, the refusal says where to put the new file to
  replace the chunk, instead of "Rename one of them".
- `index build` warns when a merge leaves out a chunk that is in the index but no
  longer complete on disk (a rebuild killed in the middle of its swap), and says so
  when it finishes a substructure build that an earlier run started.
- `index build --sdf` skips records whose ID the index cannot hold (#105). FPSim2 keeps
  IDs as 32-bit numbers, so an in-house registry number above 4,294,967,295 was indexed
  under that number modulo 2^32: a search returned another compound's ID, with its
  PubChem link, at similarity 1.0. Such records are now counted as skipped, also when
  the title line holds another number. IDs must be plain digits: `²` and other
  non-ASCII digits are no longer read as numbers, here or in `pubchem://compound/{cid}`;
  leading zeros are fine. Upgrading does not repair an index built from such a file:
  `index build --sdf` on the unchanged file keeps its chunk. Save the file again (or
  `touch` it) so the chunk is rebuilt, or build in a new `--home`.
- `index build --sdf` says why records were skipped when their IDs are the reason,
  naming the first one. ChEMBL, ZINC and vendor files, whose IDs are not numbers, used
  to end in "0 compounds" and "No built chunks" with nothing else.
- A name that is not UTF-8 on a record's title line (a Windows-1252 `±` or `µ`) no
  longer stops `index build --sdf` when the ID is in the `PUBCHEM_COMPOUND_CID` tag;
  a field that cannot be read is left empty.
- A cut-off or damaged `.sdf.gz` is reported as that, with the file's name, instead of
  as `<Boost.Python.function object at 0x...> returned a result with an exception set`.
- Substructure search is no longer dropped when its merge is interrupted (#120). The
  index goes without it while the pattern database is merged; a build stopped there
  used to come back without substructure search, call itself up to date, and delete
  the per-chunk pattern data on the next build that added compounds. Running the build
  again now finishes it, and only `--without-substructure` removes pattern data
  (including what an unfinished `--with-substructure` left behind). A failed pattern
  merge no longer says the index is unchanged.
- `pubchem_xrefs` keeps the registry IDs and PubMed articles when the patent request
  fails for any reason, not only a timeout: a dropped connection, a 500 or 502, or a
  cut-off answer used to fail the whole call. The note says why patents are missing,
  and that for the most-patented compounds (water, ethanol) PubChem's 30-second limit
  is reached every time.
- Two messages told people to run a bare `pubchem-finder ...` command, which most
  users, who start the server through `uvx`, do not have: the note on searches of the
  demo set, and the advice after an RDKit or FPSim2 upgrade. Both now give the full
  `uvx` command for the running version, the note names the command for the search
  that ran (similarity or substructure), and the upgrade advice also fits indexes
  built from your own files.
- The Caddy example in docs/clients.md could not work: with the server on its default
  loopback bind, Caddy passed the public `Host` header through and every MCP request
  was refused with a 421, while the suggested `/health` check still passed. The example
  now rewrites `Host`, and the check is an MCP request sent through the proxy with
  `curl`. The page also says how clients send the proxy's credentials and that clients
  running in a web page are not supported through it.

### Changed

- README: a "What leaves your machine" section says which tools send what to PubChem
  and which send nothing, and which of the built-in prompts point the assistant at
  PubChem. docs/clients.md marks sharing one server over HTTP as experimental and
  points to #125. ROADMAP.md follows the new milestones.
- Requests to PubChem no longer return NCBI's session cookie. A long-running server
  shares one HTTP client, so every later request of that server, for every user of a
  shared one, went out under the same session ID.
- The `find_analogues` prompt offers PubChem's live search when only the demo set was
  searched and asks first, instead of running it: that search sends the structure to
  PubChem.
- docs/clients.md and `examples/mcp_client.py` build the first index with
  `--with-substructure`, as the README does; without it, substructure search stopped
  working once the index replaced the demo set.

## [2.2.1] - 2026-11-04

### Fixed

- Parallel local searches no longer risk corrupting memory. FPSim2 releases its reference
  to an in-memory fingerprint array without holding the GIL, so two threads searching the
  same index at once (parallel MCP tool calls, or several clients over HTTP) raced on the
  array's reference count, the same fault as chembl/FPSim2#140. In-memory similarity
  searches and substructure screens now run one at a time per index.
- Oversized structures are refused with an error that names the limit, instead of
  crashing or stalling the server: SMILES and SMARTS are capped at 10,000 characters and
  1,000 atoms (the largest drug in the demo set has 492). An 80,000-character chain
  overflowed RDKit's recursive SMILES writer and killed the MCP server process. A
  similarity-map request draws at most 2,500 atoms, and molecules over 100 atoms keep
  RDKit's own layout instead of being aligned to the shared core: a 600-atom chain took
  41 s with the GIL held, which stalled every client of an HTTP server. Scaffold grouping
  skips polymer-sized hits.
- Batch mode no longer destroys its input. `props -i compounds.csv -o compounds.csv`
  (and `standardize`, `triage`, `search`, `substructure`) emptied the file before reading
  it. Writing over the input is now refused, and the input is opened and its header
  checked before the output file is touched, so a bad input leaves an existing output
  file as it was.
- `index build --sdf` no longer loses compounds when files share a name. Each file becomes
  a chunk named after it, so `--sdf vendor_a/library.sdf --sdf vendor_b/library.sdf`
  kept only one of them, while reporting both as built, and a file with the same name
  in a later run was skipped as already built. Both are now refused with a message that
  says to rename one. Chunks record the file they were built from, and running the
  command again after a file changed rebuilds its chunk instead of skipping it.
- `pubchem_xrefs` no longer fails for well-studied compounds when PubChem runs out of
  time listing their patents (#90): patents, most of the response (3.1 of 3.6 MB for
  aspirin), are fetched on their own, and if that request times out the registry IDs and
  PubMed articles are still returned, with `unavailable: ["PatentID"]` and a `note`.
- A PubChem timeout on a lookup suggests a retry; only similarity and substructure
  searches suggest a narrower query. Timeout messages no longer call every request a
  search.

## [2.2.0] - 2026-11-02

### Added

- Local similarity and substructure search work on first use: until an index is built,
  they search a built-in demo set of 1,932 FDA-approved drugs (PubChem's parent compounds
  of what it links to Drugs@FDA, public domain), built into the data directory in about a
  second on the first search. Results say so with `search_scope: "demo"` and a `note`
  that points to PubChem's live searches and the build command; `index_status` reports
  `search_scope` and `demo_compounds`. Once an index is built, searches use it without a
  restart. `index build --demo` builds the demo index ahead of time and
  `PUBCHEM_FINDER_DEMO=0` turns the fallback off. `scripts/build_demo_set.py` regenerates
  the data.
- Three MCP prompts start common workflows by name (in Claude Code,
  `/pubchem-finder:compound_dossier aspirin`): `compound_dossier` (a cited profile:
  identity, structure, safety, pharmacology, bioassay targets), `find_analogues` (close
  analogues grouped by scaffold and drawn with what changed, optional `threshold`) and
  `triage_structures` (medicinal-chemistry triage of comma-separated SMILES, ranked).
- Two MCP resources: `pubchem://compound/{cid}` (a PubChem compound record as JSON; in
  Claude Code, `@pubchem-finder:pubchem://compound/2244`) and `pubchem://index/status`.
  Errors say what went wrong: an unknown CID is "not found", other failures name their
  cause.
- The index manifest records an index format version. A pubchem-finder that meets an
  index written in a newer format says so and what to do, instead of misreading it.
  Indexes from 2.1 and earlier are format 1 and keep working.

### Changed

- When there is no local index, the message (and `index_status`) gives a command that
  runs from any terminal: `uvx --python 3.14 pubchem-finder==<version> index build --limit 1
  --with-substructure`, with `--home` when the server uses another data directory, and
  what it costs (PubChem's first chunk is now a 323 MB download). Plugin and MCP users
  have no `pubchem-finder` command on their PATH, so the old hint didn't run for them.
  `index build` shows download progress.
- `similarity_search`, `substructure_search` and `pubchem_annotations` ask Claude Code to
  show results of up to 100,000 characters inline (`anthropic/maxResultSizeChars`); by
  default it saves results over 50,000 characters to a file instead, which large searches
  can pass.
- The Claude plugin runs 2.1.1 and passes `--python=3.14` to uvx, like the other
  documented launch commands, so plugin users get the tool annotations and stay off
  Python 3.15.

### Fixed

- In-memory searches no longer risk failing with "Chunk N worker died: ...
  GenericSearch(): incompatible function arguments" or crashing the server. FPSim2 0.7.4
  updates the fingerprint arrays' reference counts unsafely when it splits a search
  across threads, and can free them mid-search; searches now run in the calling thread.
  Seen on Python 3.11, where it made long sessions and batch searches fail; reference
  counts drift on 3.12 too.
- The declared minimum versions work again: Pillow 10.1 (`draw_molecules` needs
  `ImageFont.load_default(size=...)`; with Pillow 10.0 it failed) and typer 0.16 (0.15
  crashed on usage errors with click 8.2 or later). CI now runs the tests against the
  lowest allowed version of every direct dependency.

## [2.1.1] - 2026-10-28

### Added

- Every MCP tool has a title and annotations (read-only, non-destructive, idempotent;
  open-world for the `pubchem_*` tools, which call PubChem), as Anthropic's Software
  Directory Policy requires. The server reports its title, website and icon, and the
  MCP registry entry shows the icon.
- docs/clients.md has a tested setup for VS Code with GitHub Copilot, the first from an
  outside contributor (@ashwinruke, [#37](https://github.com/MysterionRise/pubchem-finder/pull/37)).
  Its notes cover a Copilot quirk worth knowing: the model can show a `similarity_search`
  step and still report compounds and scores the tool never returned.

### Changed

- The documented launch commands, the example client and the MCP registry entry pass
  `--python 3.14` to uvx. uvx uses uv's newest Python whatever the package supports, and
  on Python 3.15 FPSim2, numexpr and ndindex publish no wheels yet, so installing
  compiles them (about half a minute with a C/C++ compiler) or fails without one. CI
  checks that dependencies install from wheels on Python 3.11 and 3.14 on Linux, macOS
  and Windows, and a weekly job reports when Python 3.15 catches up.
- The Claude plugin runs one pinned release (`uvx pubchem-finder==2.1.0 mcp`), ships a
  lockfile (`plugins/pubchem-finder/pyproject.toml` and `uv.lock`: the exact dependency
  set, from PyPI with hashes) and has an icon and documentation and support links, as
  Anthropic's plugin directory requires. `make plugin-lock VERSION=...` points it at a new
  release; tests keep its files in step and never ahead of the package.

### Fixed

- Extending an index after an RDKit or FPSim2 upgrade no longer leaves it unusable:
  `index build` checks which versions fingerprinted the existing chunks before
  downloading anything and says what to do, and a failed merge keeps the previous index
  searchable (the merged file replaces it only when complete). Chunks record their RDKit
  and FPSim2 versions, and pattern chunks are rebuilt when FPSim2 changes.
- Indexes whose chunks overlap (for example an `--sdf` file repeating compounds of a
  PubChem chunk) no longer return the same CID twice from `similarity_search` or
  `substructure_search`.
- Unexpected failures in MCP tools (a damaged cache or index file, a disk error) reach
  the assistant as a short error naming the problem instead of a bare "Error
  executing tool". The PubChem response cache falls back to the network when it can't
  be read or written, and caches only JSON objects. `index status` and `index build`
  report such errors without a traceback.

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

[Unreleased]: https://github.com/MysterionRise/pubchem-finder/compare/v2.2.3...HEAD
[2.2.3]: https://github.com/MysterionRise/pubchem-finder/releases/tag/v2.2.3
[2.2.2]: https://github.com/MysterionRise/pubchem-finder/releases/tag/v2.2.2
[2.2.1]: https://github.com/MysterionRise/pubchem-finder/releases/tag/v2.2.1
[2.2.0]: https://github.com/MysterionRise/pubchem-finder/releases/tag/v2.2.0
[2.1.1]: https://github.com/MysterionRise/pubchem-finder/releases/tag/v2.1.1
[2.1.0]: https://github.com/MysterionRise/pubchem-finder/releases/tag/v2.1.0
[2.0.1]: https://github.com/MysterionRise/pubchem-finder/releases/tag/v2.0.1
[2.0.0]: https://github.com/MysterionRise/pubchem-finder/releases/tag/v2.0.0
