# Releasing

How a version goes out: PyPI, the tag and GitHub Release, the MCP registry, then the
Claude plugin. `.github/workflows/release.yml` does the publishing; this page is the
order of steps and what to do when one fails.

## Steps

1. **Prepare the files.** On the branch at the top of the release's stack (or a new
   `claude/release-x.y.z` branch):

   ```bash
   make release-prep VERSION=x.y.z
   ```

   It moves the CHANGELOG's Unreleased entries under `## [x.y.z] - <today>`, fixes the
   compare links, sets the version in `pyproject.toml`, `server.json` and
   `CITATION.cff` (with `date-released`) and runs `uv lock`. Review the diff, commit,
   open the PR. If the release slips to another day, update the two dates.

2. **Dry run** from that branch. It builds, tests and checks the package, and builds a
   real index from PubChem's first SDF chunk (about 10 minutes in all):

   ```bash
   gh workflow run release.yml --ref claude/release-x.y.z -f version=x.y.z
   ```

3. **Merge** the stack: `gh stack merge <stack>` (or the single PR).

4. **Release** from master:

   ```bash
   gh workflow run release.yml -f version=x.y.z -f dry_run=false
   ```

   The run waits for CI on master's commit, tags it, drafts the GitHub Release, then
   stops at the `pypi` environment: approve it under the run's **Review deployments**.
   It then publishes to PyPI, publishes the GitHub Release and lists the version in the
   MCP registry.

5. **Check** what users get:

   ```bash
   uvx --python 3.14 pubchem-finder@x.y.z --version
   curl -s "https://registry.modelcontextprotocol.io/v0/servers?search=pubchem-finder" | jq '.servers[].server.version'
   ```

6. **Plugin.** The run's `plugin-branch` job pushes `release/plugin-x.y.z` with the
   plugin pointed at the release and a fresh `uv.lock`; its summary has the link to open
   the PR. Merging that PR fast-forwards `plugin-release`, the branch Anthropic's plugin
   directory follows (`promote-plugin.yml`), and the directory scans it as a new version.
   Each version waits for a directory reviewer (the plugin launches a package with uvx),
   so promote every minor release, and a patch only when plugin users gain something:
   a fix in a tool, new annotations, a security fix. An unpromoted plugin stays on its
   older pin, which the tests allow. By hand instead: `make plugin-lock VERSION=x.y.z`.

## When a step fails

| What failed | What to do |
|---|---|
| Guard: not on master | Real releases only run from master; merge first. |
| Guard: no successful CI | Start it with `gh workflow run ci.yml --ref master` and run the release again once it passes. |
| Versions, changelog and citation | `make release-prep` again, or fix what the error names. |
| Lint, tests, package metadata, dry-run end-to-end | Fix, push a new commit, start again. |
| Tag already points elsewhere | Only if PyPI doesn't have the version: delete the tag (`gh api -X DELETE repos/MysterionRise/pubchem-finder/git/refs/tags/vx.y.z`; admins only) and start again. |
| `publish`, `github-release`, `mcp-registry` or `plugin-branch` | **Re-run failed jobs** on the same run: every job resumes safely (`skip-existing` on PyPI, the tag and draft are reused). |
| PyPI approval rejected | Re-run the job to approve later, or delete the draft release and the tag. |
| The plugin PR's lock is over 256 KiB | Restrict `[tool.uv] environments` in `plugins/pubchem-finder/pyproject.toml` (see `tests/unit/test_metadata.py`). |
