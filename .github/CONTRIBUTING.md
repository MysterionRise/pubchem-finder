# Contributing

Thanks for helping with PubChem Finder. The project is small; this page is the whole process.

## Setup

Install [uv](https://docs.astral.sh/uv/), then:

```bash
git clone https://github.com/MysterionRise/pubchem-finder.git
cd pubchem-finder
uv sync          # creates .venv with the project and dev tools (RDKit comes from PyPI wheels)
make check       # ruff lint + format check, mypy --strict, pytest with coverage
```

Python 3.11 to 3.14 are supported and tested in CI.

## Workflow

1. Branch from `master`.
2. Make the change with tests. Unit tests go in `tests/unit/`, end-to-end flows in
   `tests/integration/`, and tests that start the real MCP server process in `tests/e2e/`.
   Tests must run offline and finish in seconds; mark anything that calls PubChem with
   `@pytest.mark.live`.
3. Run `make check`. CI runs the same commands, so a green local run means a green PR.
4. Open a pull request. Keep it focused and describe what changed and why.

## Good first contributions

- **Set it up in your MCP client and write down what worked.** The server runs in any
  client that launches local MCP servers, but only Claude's setup is documented.
  [docs/clients.md](../docs/clients.md) has the settings, four questions that check the
  setup, and the format for a new section; no Python changes needed.
- Issues labelled
  [good first issue](https://github.com/MysterionRise/pubchem-finder/labels/good%20first%20issue)
  are scoped to be done in an evening.

## What CI checks

- **CI** (every pull request, whatever its base, and every push to master): lint, format
  and `mypy --strict`; the test suite on Linux (Python 3.11 to 3.14), macOS and Windows;
  a `lowest` job on the oldest version of every dependency that `pyproject.toml` allows;
  a `wheels` job that checks every dependency installs from wheels on Python 3.11 and
  3.14; and a `package` job on all three systems that builds the wheel, installs it into
  a clean environment with the newest dependencies from PyPI, and runs `scripts/smoke.py`
  against the installed command. `ci-ok` passes when all of them do; it is the check
  master requires.
- **Next Python** (weekly): whether the latest release installs from wheels on Python
  3.15 yet, and which packages still lack them.
- **Live PubChem** (weekly and on demand): the `live` tests against the real PubChem
  services. Run them locally with `PUBCHEM_LIVE=1 uv run pytest -m live`.
- **Monthly real-data check**: the real-data end-to-end workflow below, on master. A
  failing scheduled run (this one or the weekly live tests) opens an issue labelled
  `ci-failure`, which closes itself when the next run passes.
- **Real-data end-to-end** (on demand, monthly, and in every release dry run): downloads real PubChem SDF
  chunks, builds an index with substructure search, runs every CLI command and every MCP
  tool against it and against PubChem, and writes build time, file sizes, memory and query
  timings to the job summary. Locally: `uv run pubchem-finder index build --limit 1
  --with-substructure`, then `uv run python scripts/smoke.py --home <data dir> --live`.

## Releasing

[docs/releasing.md](../docs/releasing.md) has the steps and what to do when one
fails. In short: `make release-prep VERSION=x.y.z`, a dry run, merge, then the release
workflow from master, then the plugin PR it prepares.

## Code style

- `ruff format` decides formatting; `ruff check` decides lint. Both run in CI and as
  pre-commit hooks (`uv run pre-commit install`), with zizmor and actionlint for the
  workflows; `make check` and CI run them all.
- `mypy --strict` must pass. RDKit is untyped; keep its objects behind small, typed
  wrappers rather than sprinkling `Any` through the codebase.
- Public functions get a short docstring. Comment the why, not the what.

## Dependencies

Add runtime dependencies with `uv add <package>` and dev tools with `uv add --dev <package>`,
and commit the updated `uv.lock`. Do not add a dependency the code does not import yet.

## Reporting problems

Open an issue with the SMILES or input that triggers it, the command you ran, and the
versions printed by `uv run python -c "import rdkit, pubchem_finder; print(rdkit.__version__, pubchem_finder.__version__)"`.
