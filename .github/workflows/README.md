# GitHub Actions CI/CD Pipeline

This directory contains the automated CI/CD workflows for PubChem Finder.

---

## 📋 CI Workflow (`ci.yml`)

The main CI pipeline runs on:
- **Push** to `main`, `develop`, or any `feature/*` branch
- **Pull requests** to `main` or `develop`

### Pipeline Stages

#### 1. **Test Job** (Matrix: Python 3.10, 3.11)

Runs on: `ubuntu-latest`

**Steps:**

1. ✅ **Setup**
   - Checkout code
   - Set up Python (3.10, 3.11)
   - Install Poetry
   - Cache dependencies

2. ✅ **Validation**
   - Run `scripts/validate.py`
   - Checks: directory structure, Python syntax, imports, test structure

3. ✅ **Code Quality**
   - **Black**: Code formatting check
   - **Ruff**: Linting
   - **Mypy**: Type checking (non-blocking)

4. ✅ **Unit Tests**
   - Run `tests/unit/` with coverage
   - 50+ tests covering core modules

5. ✅ **Integration Tests**
   - Run `tests/integration/` with coverage
   - 15+ end-to-end workflow tests

6. ✅ **Performance Tests** (Non-blocking)
   - Run `tests/test_performance.py`
   - Benchmarks with timing info

7. ✅ **Coverage Report**
   - Generate final coverage (XML + term)
   - Upload HTML report as artifact

#### 2. **Lint Job** (Python 3.11)

Separate linting job:
- Black formatting check
- Ruff linting

---

## 🎯 What Gets Tested

### Code Quality Checks
```yaml
✓ Code formatting (black --check)
✓ Linting (ruff check)
✓ Type hints (mypy)
```

### Test Suite (70+ tests)
```yaml
✓ Unit tests (50+)
✓ Integration tests (15+)
✓ Performance benchmarks (5+)
```

### Validation
```yaml
✓ Project structure
✓ Python syntax
✓ Import integrity
✓ Test structure
```

---

## 📊 CI Execution Flow

```
┌─────────────────────────────────────────┐
│  Push/PR to main/develop/feature/*      │
└────────────────┬────────────────────────┘
                 │
                 ▼
┌─────────────────────────────────────────┐
│  Matrix: Python 3.10 & 3.11             │
│  ┌───────────────────────────────────┐  │
│  │ 1. Setup & Install Dependencies   │  │
│  └───────────────────────────────────┘  │
│  ┌───────────────────────────────────┐  │
│  │ 2. Validate Project Structure     │  │
│  │    python scripts/validate.py     │  │
│  └───────────────────────────────────┘  │
│  ┌───────────────────────────────────┐  │
│  │ 3. Code Quality Checks            │  │
│  │    - black --check                │  │
│  │    - ruff check                   │  │
│  │    - mypy src                     │  │
│  └───────────────────────────────────┘  │
│  ┌───────────────────────────────────┐  │
│  │ 4. Run Unit Tests                 │  │
│  │    pytest tests/unit --cov        │  │
│  └───────────────────────────────────┘  │
│  ┌───────────────────────────────────┐  │
│  │ 5. Run Integration Tests          │  │
│  │    pytest tests/integration       │  │
│  └───────────────────────────────────┘  │
│  ┌───────────────────────────────────┐  │
│  │ 6. Run Performance Tests          │  │
│  │    pytest test_performance.py     │  │
│  └───────────────────────────────────┘  │
│  ┌───────────────────────────────────┐  │
│  │ 7. Generate Coverage Report       │  │
│  │    Upload artifacts               │  │
│  └───────────────────────────────────┘  │
└─────────────────────────────────────────┘
```

---

## ✅ Success Criteria

For CI to pass, all of the following must succeed:

1. ✅ **Validation**: Project structure is valid
2. ✅ **Formatting**: Code follows Black style
3. ✅ **Linting**: No Ruff errors
4. ✅ **Unit Tests**: All 50+ tests pass
5. ✅ **Integration Tests**: All 15+ tests pass
6. ✅ **Coverage**: Maintained at 85%+

**Non-blocking** (informational only):
- ⚠️ Type checking (mypy)
- ⚠️ Performance tests

---

## 📦 Artifacts

After each CI run, the following artifacts are available:

- **Coverage Report** (HTML)
  - Path: `htmlcov/`
  - Available for: 90 days
  - Per Python version (3.10, 3.11)

---

## 🔍 Viewing CI Results

### On GitHub

1. Go to repository → **Actions** tab
2. Select a workflow run
3. View job logs and test results
4. Download coverage artifacts

### Locally (Before Pushing)

Run the same checks that CI runs:

```bash
# All quality checks
make quality

# All tests
make test

# Validate structure
python scripts/validate.py

# Everything
make quality && make test && python scripts/validate.py
```

---

## 🐛 Debugging CI Failures

### Common Issues

**1. Import Errors**
```
Solution: Check pyproject.toml dependencies
Verify: python scripts/validate.py
```

**2. Test Failures**
```
Solution: Run tests locally first
Command: make test
```

**3. Formatting Issues**
```
Solution: Auto-format code
Command: make format
```

**4. Linting Errors**
```
Solution: Fix with ruff
Command: poetry run ruff check --fix src tests
```

---

## 🚀 Local Testing Matches CI

The CI environment uses:
- Ubuntu latest
- Python 3.10 & 3.11
- Poetry 1.7.1
- RDKit via conda

To match locally:
```bash
# Use Docker (exact same env as CI)
make docker-test

# Or install dependencies
poetry install
make test
```

---

## 📈 CI Performance

Typical CI run times:

| Stage | Duration |
|-------|----------|
| Setup & Install | ~2-3 min |
| Validation | ~5 sec |
| Code Quality | ~10 sec |
| Unit Tests | ~30 sec |
| Integration Tests | ~20 sec |
| Performance Tests | ~15 sec |
| Coverage Report | ~5 sec |
| **Total** | **~4-5 min** |

---

## 🔧 Modifying CI

### Adding New Test Categories

Edit `.github/workflows/ci.yml`:

```yaml
- name: Run new test category
  run: poetry run pytest tests/new_category -v
```

### Changing Python Versions

Edit the matrix:

```yaml
strategy:
  matrix:
    python-version: ["3.10", "3.11", "3.12"]  # Add 3.12
```

### Adding New Quality Checks

```yaml
- name: Run security check
  run: poetry run bandit -r src
```

---

## 📝 CI Badge

Add to README.md:

```markdown
[![CI](https://github.com/MysterionRise/pubchem-finder/workflows/CI/badge.svg)](https://github.com/MysterionRise/pubchem-finder/actions)
```

---

## 🎓 Best Practices

1. **Run tests locally first**
   ```bash
   make quality && make test
   ```

2. **Check CI before merging**
   - Wait for green checkmark
   - Review any warnings

3. **Keep CI fast**
   - Use caching
   - Parallel jobs where possible
   - Non-blocking for non-critical checks

4. **Monitor coverage**
   - Download artifacts
   - Check coverage trends
   - Maintain 85%+ threshold

---

## 🔗 Related Files

- `pyproject.toml` - Dependencies and tool config
- `Makefile` - Local commands that mirror CI
- `scripts/validate.py` - Pre-test validation
- `TESTING.md` - Complete testing guide

---

**CI Status**: ✅ Fully Configured | 70+ Tests | Multi-Python | Coverage Tracking
