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

Python 3.11, 3.12 and 3.13 are supported and tested in CI.

## Workflow

1. Branch from `master`.
2. Make the change with tests. Unit tests go in `tests/unit/`, end-to-end flows in
   `tests/integration/`. Tests must run offline; mark anything slow with `@pytest.mark.slow`.
3. Run `make check`. CI runs the same commands, so a green local run means a green PR.
4. Open a pull request. Keep it focused and describe what changed and why.

## Code style

- `ruff format` decides formatting; `ruff check` decides lint. Both run in CI and as
  pre-commit hooks (`uv tool install pre-commit && uv run pre-commit install`).
- `mypy --strict` must pass. RDKit is untyped; keep its objects behind small, typed
  wrappers rather than sprinkling `Any` through the codebase.
- Public functions get a short docstring. Comment the why, not the what.

## Dependencies

Add runtime dependencies with `uv add <package>` and dev tools with `uv add --dev <package>`,
and commit the updated `uv.lock`. Do not add a dependency the code does not import yet.

## Reporting problems

Open an issue with the SMILES or input that triggers it, the command you ran, and the
versions printed by `uv run python -c "import rdkit, pubchem_finder; print(rdkit.__version__, pubchem_finder.__version__)"`.
