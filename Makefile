.PHONY: help install lint format typecheck test check release-prep plugin-lock clean

help:  ## Show available targets
	@awk 'BEGIN {FS = ":.*?## "} /^[a-zA-Z_-]+:.*?## / {printf "  \033[36m%-12s\033[0m %s\n", $$1, $$2}' $(MAKEFILE_LIST)

install:  ## Create the virtualenv and install the project with dev tools
	uv sync

lint:  ## Lint with ruff
	uv run ruff check .

format:  ## Format with ruff
	uv run ruff format .

typecheck:  ## Type-check with mypy
	uv run mypy

test:  ## Run the test suite with coverage
	uv run pytest --cov --cov-report=term-missing

check: lint typecheck test  ## Everything CI runs
	uv run ruff format --check .

release-prep:  ## Roll the changelog and set every version: make release-prep VERSION=2.1.2
	uv run --no-project python scripts/release_prep.py prep $(VERSION)

plugin-lock:  ## Point the Claude plugin at a release (the release workflow does it): make plugin-lock VERSION=x.y.z
	uv run python scripts/plugin_release.py $(VERSION)

clean:  ## Remove caches and build artifacts
	rm -rf .pytest_cache .mypy_cache .ruff_cache htmlcov coverage.xml .coverage dist build
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
