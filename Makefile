.PHONY: help install test test-watch lint format type-check clean docker-build docker-test

help:  ## Show this help message
	@echo 'Usage: make [target]'
	@echo ''
	@echo 'Available targets:'
	@awk 'BEGIN {FS = ":.*?## "} /^[a-zA-Z_-]+:.*?## / {printf "  \033[36m%-20s\033[0m %s\n", $$1, $$2}' $(MAKEFILE_LIST)

install:  ## Install dependencies with Poetry
	poetry install

install-dev:  ## Install all dependencies including dev tools
	poetry install --with dev

test:  ## Run tests with coverage
	poetry run pytest --cov --cov-report=term-missing --cov-report=html -v

test-unit:  ## Run only unit tests
	poetry run pytest tests/unit -v

test-integration:  ## Run only integration tests
	poetry run pytest tests/integration -v

test-watch:  ## Run tests in watch mode
	poetry run pytest-watch -- --cov -v

test-performance:  ## Run performance benchmarks
	poetry run pytest tests/test_performance.py -v --durations=10

lint:  ## Run linting with ruff
	poetry run ruff check src tests

format:  ## Format code with black
	poetry run black src tests

format-check:  ## Check code formatting
	poetry run black --check src tests

type-check:  ## Run type checking with mypy
	poetry run mypy src

pre-commit:  ## Run pre-commit hooks on all files
	poetry run pre-commit run --all-files

quality:  ## Run all quality checks
	@echo "Running code formatting..."
	poetry run black src tests
	@echo "Running linter..."
	poetry run ruff check src tests
	@echo "Running type checker..."
	poetry run mypy src
	@echo "All quality checks passed!"

clean:  ## Clean up cache and build files
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name "*.pyc" -delete
	find . -type f -name "*.pyo" -delete
	find . -type d -name "*.egg-info" -exec rm -rf {} + 2>/dev/null || true
	find . -type d -name ".pytest_cache" -exec rm -rf {} + 2>/dev/null || true
	find . -type d -name ".mypy_cache" -exec rm -rf {} + 2>/dev/null || true
	find . -type d -name ".ruff_cache" -exec rm -rf {} + 2>/dev/null || true
	rm -rf htmlcov/
	rm -rf dist/
	rm -rf build/
	rm -f .coverage

docker-build:  ## Build Docker image
	docker-compose build

docker-test:  ## Run tests in Docker
	docker-compose run --rm test

docker-shell:  ## Open shell in Docker container
	docker-compose run --rm dev

docker-clean:  ## Clean Docker images and containers
	docker-compose down -v
	docker rmi pubchem-finder:dev 2>/dev/null || true

coverage:  ## Generate and open coverage report
	poetry run pytest --cov --cov-report=html
	@echo "Opening coverage report..."
	@python -m webbrowser htmlcov/index.html 2>/dev/null || echo "Please open htmlcov/index.html manually"

init:  ## Initialize project (install dependencies and pre-commit hooks)
	poetry install
	poetry run pre-commit install
	@echo "Project initialized! Run 'make test' to verify everything works."

validate:  ## Validate project structure and imports
	@echo "Validating project structure..."
	@python scripts/validate.py

benchmark:  ## Run performance benchmarks and show results
	poetry run pytest tests/test_performance.py -v --durations=10 --tb=short
