# Contributing to PubChem Finder

Thank you for considering contributing to PubChem Finder! This guide will help you get started.

---

## 🚀 Quick Start

### 1. Fork & Clone

```bash
git clone https://github.com/YOUR_USERNAME/pubchem-finder.git
cd pubchem-finder
```

### 2. Set Up Development Environment

**Option A: Docker (Recommended)**
```bash
make docker-build
make docker-shell
```

**Option B: Local with Poetry**
```bash
poetry install
poetry shell
```

### 3. Create a Branch

```bash
git checkout -b feature/my-awesome-feature
```

### 4. Make Changes

Edit code, write tests, update docs.

### 5. Run Quality Checks

```bash
# All checks
make quality

# All tests
make test

# Or use Docker
make docker-test
```

### 6. Commit & Push

```bash
git add .
git commit -m "Add awesome feature"
git push origin feature/my-awesome-feature
```

### 7. Create Pull Request

Open a PR on GitHub with a clear description.

---

## 📋 Development Workflow

### Before You Start

1. Check existing issues and PRs
2. Discuss major changes in an issue first
3. Read the implementation plan (`IMPLEMENTATION_PLAN.md`)

### Code Standards

- **Style**: Black formatter (line length: 88)
- **Linting**: Ruff with standard rules
- **Type Hints**: 100% coverage required
- **Docstrings**: All public APIs
- **Tests**: Required for all new features

### Running Checks

```bash
# Format code
make format

# Lint
make lint

# Type check
make type-check

# All quality checks
make quality
```

---

## 🧪 Testing Requirements

### Writing Tests

All new code must include tests:

1. **Unit tests** for individual functions
2. **Integration tests** for workflows
3. **Performance tests** if relevant

### Test Structure

```python
# tests/unit/test_my_feature.py
import pytest
from pubchem_finder import MyFeature


class TestMyFeature:
    """Tests for MyFeature."""

    def test_basic_functionality(self):
        """Test basic use case."""
        feature = MyFeature()
        result = feature.do_something()
        assert result is not None

    def test_error_handling(self):
        """Test error cases."""
        with pytest.raises(ValueError):
            MyFeature(invalid_param=True)
```

### Running Tests

```bash
# All tests
make test

# Specific file
poetry run pytest tests/unit/test_my_feature.py -v

# With coverage
poetry run pytest --cov -v
```

### Coverage Requirement

- Maintain **85%+** overall coverage
- New code should have **>90%** coverage

---

## 📝 Documentation

### Update Documentation When:

- Adding new features → Update README.md
- Changing APIs → Update docstrings
- Adding tests → Update TESTING.md
- Changing workflow → Update CONTRIBUTING.md

### Docstring Format

```python
def calculate_similarity(fp1: np.ndarray, fp2: np.ndarray) -> float:
    """
    Calculate Tanimoto similarity between two fingerprints.

    Args:
        fp1: First fingerprint (binary array)
        fp2: Second fingerprint (binary array)

    Returns:
        Similarity score between 0 and 1

    Raises:
        ValueError: If fingerprints have different lengths

    Example:
        >>> fp1 = np.array([1, 0, 1])
        >>> fp2 = np.array([1, 1, 0])
        >>> calculate_similarity(fp1, fp2)
        0.333
    """
```

---

## 🔄 Pull Request Process

### Before Submitting

1. ✅ Run all quality checks: `make quality`
2. ✅ Run all tests: `make test`
3. ✅ Update documentation
4. ✅ Add yourself to contributors (if first PR)

### PR Description Template

```markdown
## Description
Brief description of changes

## Type of Change
- [ ] Bug fix
- [ ] New feature
- [ ] Breaking change
- [ ] Documentation update

## Testing
- [ ] Unit tests added/updated
- [ ] Integration tests added/updated
- [ ] All tests pass locally
- [ ] Coverage maintained/improved

## Checklist
- [ ] Code follows style guidelines
- [ ] Self-review completed
- [ ] Documentation updated
- [ ] No new warnings
```

### Review Process

1. **Automated Checks**: CI must pass
2. **Code Review**: Maintainer review
3. **Testing**: Verify all tests pass
4. **Documentation**: Check docs are updated
5. **Merge**: Squash and merge

---

## 🐛 Reporting Bugs

### Bug Report Template

```markdown
**Describe the bug**
Clear description of the issue

**To Reproduce**
Steps to reproduce:
1. Import '...'
2. Call function '...'
3. See error

**Expected behavior**
What you expected to happen

**Environment**
- OS: [e.g., Ubuntu 22.04]
- Python: [e.g., 3.11]
- PubChem Finder: [e.g., 2.0.0]

**Additional context**
Any other relevant information
```

---

## 💡 Suggesting Features

### Feature Request Template

```markdown
**Is your feature related to a problem?**
Description of the problem

**Proposed solution**
How you'd like it to work

**Alternatives considered**
Other approaches you've thought about

**Phase**
Which implementation phase does this fit? (1-5)
```

---

## 🏗️ Project Structure

```
pubchem-finder/
├── src/pubchem_finder/     # Main source code
│   ├── core/               # Core abstractions
│   ├── chemistry/          # Chemistry utilities
│   ├── fingerprints/       # Fingerprint generators
│   └── ...                 # Future modules
├── tests/                  # All tests
│   ├── unit/              # Unit tests
│   ├── integration/       # Integration tests
│   └── fixtures/          # Test data
├── scripts/               # Utility scripts
└── docs/                  # Documentation
```

---

## 🎯 Areas to Contribute

### Phase 1 (Current) - Enhancements

- [ ] Add more fingerprint types (RDKit, AtomPair)
- [ ] Improve error messages
- [ ] Add more molecular properties
- [ ] Performance optimizations
- [ ] Better test coverage

### Phase 2 (Next) - Similarity Search

- [ ] Implement search engine
- [ ] PostgreSQL integration
- [ ] Vector storage
- [ ] CLI commands

### Phase 3+ (Future)

- [ ] REST API
- [ ] Web interface
- [ ] ML-based embeddings
- [ ] Advanced search algorithms

---

## 📚 Resources

### Documentation

- [README.md](../README.md) - Project overview
- [IMPLEMENTATION_PLAN.md](../IMPLEMENTATION_PLAN.md) - Roadmap
- [TESTING.md](../TESTING.md) - Testing guide
- [REPLATFORMING_PROMPT.md](../REPLATFORMING_PROMPT.md) - Technical spec

### External Resources

- [RDKit Documentation](https://www.rdkit.org/docs/)
- [Python Type Hints](https://docs.python.org/3/library/typing.html)
- [Pytest Documentation](https://docs.pytest.org/)

---

## 🤝 Code of Conduct

### Our Standards

- Be respectful and inclusive
- Provide constructive feedback
- Focus on what's best for the project
- Show empathy towards other contributors

### Unacceptable Behavior

- Harassment or discriminatory language
- Personal attacks
- Publishing others' private information
- Other unprofessional conduct

---

## 📧 Getting Help

- **Questions**: Open a GitHub issue
- **Discussions**: Use GitHub Discussions
- **Bugs**: File a bug report
- **Security**: Email maintainers directly

---

## ✅ Contribution Checklist

Before submitting:

- [ ] Code follows style guide (Black + Ruff)
- [ ] All tests pass (`make test`)
- [ ] Coverage maintained at 85%+
- [ ] Type hints added for new code
- [ ] Docstrings added for public APIs
- [ ] Documentation updated
- [ ] No merge conflicts
- [ ] Commit messages are clear
- [ ] PR description is complete

---

## 🎉 Recognition

Contributors will be:
- Listed in CONTRIBUTORS.md
- Mentioned in release notes
- Credited in documentation

---

**Thank you for contributing to PubChem Finder!** 🚀
