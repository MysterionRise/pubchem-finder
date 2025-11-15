#!/usr/bin/env python3
"""
Validation script to check project structure and basic imports.

This script verifies:
1. All required directories exist
2. All Python files have valid syntax
3. Core modules can be imported (if dependencies available)
4. Test files are properly structured
"""

import ast
import sys
from pathlib import Path
from typing import List, Tuple


def check_directory_structure() -> List[str]:
    """Verify all required directories exist."""
    errors = []
    required_dirs = [
        "src/pubchem_finder",
        "src/pubchem_finder/core",
        "src/pubchem_finder/chemistry",
        "src/pubchem_finder/fingerprints",
        "tests/unit",
        "tests/integration",
        "tests/fixtures",
    ]

    for dir_path in required_dirs:
        if not Path(dir_path).is_dir():
            errors.append(f"Missing required directory: {dir_path}")

    return errors


def check_python_syntax() -> List[str]:
    """Check all Python files for syntax errors."""
    errors = []
    python_files = list(Path("src").rglob("*.py")) + list(Path("tests").rglob("*.py"))

    for file_path in python_files:
        try:
            with open(file_path, "r", encoding="utf-8") as f:
                source = f.read()
                ast.parse(source, filename=str(file_path))
        except SyntaxError as e:
            errors.append(f"Syntax error in {file_path}: {e}")
        except Exception as e:
            errors.append(f"Error reading {file_path}: {e}")

    return errors


def check_imports() -> List[Tuple[str, bool, str]]:
    """Attempt to import core modules."""
    results = []

    # Add src to path temporarily for validation
    src_path = Path(__file__).parent.parent / "src"
    if src_path.exists():
        sys.path.insert(0, str(src_path))

    # Core imports to test
    imports_to_test = [
        ("pubchem_finder", "Main package"),
        ("pubchem_finder.core", "Core module"),
        ("pubchem_finder.core.molecule", "Molecule class"),
        ("pubchem_finder.core.fingerprint", "Fingerprint base classes"),
        ("pubchem_finder.chemistry", "Chemistry module"),
        ("pubchem_finder.fingerprints", "Fingerprints module"),
    ]

    for module_name, description in imports_to_test:
        try:
            __import__(module_name)
            results.append((module_name, True, "OK"))
        except ImportError as e:
            # Expected if RDKit not installed
            if "rdkit" in str(e).lower() or "numpy" in str(e).lower():
                results.append(
                    (module_name, False, "Missing dependencies (RDKit/numpy) - expected in dev env")
                )
            else:
                results.append((module_name, False, f"Import error: {e}"))
        except Exception as e:
            results.append((module_name, False, f"Unexpected error: {e}"))

    return results


def check_test_structure() -> List[str]:
    """Verify test files follow expected structure."""
    errors = []

    # Check for conftest.py
    if not Path("tests/conftest.py").exists():
        errors.append("Missing tests/conftest.py")

    # Check unit tests
    unit_tests = list(Path("tests/unit").glob("test_*.py"))
    if not unit_tests:
        errors.append("No unit test files found in tests/unit/")

    # Check that test files have proper structure
    for test_file in unit_tests:
        try:
            with open(test_file, "r") as f:
                source = f.read()
                tree = ast.parse(source)

                # Check for test functions
                test_functions = [
                    node.name
                    for node in ast.walk(tree)
                    if isinstance(node, ast.FunctionDef) and node.name.startswith("test_")
                ]

                if not test_functions:
                    errors.append(f"No test functions found in {test_file.name}")

        except Exception as e:
            errors.append(f"Error analyzing {test_file}: {e}")

    return errors


def check_required_files() -> List[str]:
    """Check for required project files."""
    errors = []
    required_files = [
        "pyproject.toml",
        "README.md",
        ".gitignore",
        ".github/workflows/ci.yml",
        "src/pubchem_finder/__init__.py",
    ]

    for file_path in required_files:
        if not Path(file_path).exists():
            errors.append(f"Missing required file: {file_path}")

    return errors


def main() -> int:
    """Run all validation checks."""
    print("=" * 70)
    print("PubChem Finder - Project Validation")
    print("=" * 70)

    all_passed = True

    # Check directory structure
    print("\n1. Checking directory structure...")
    dir_errors = check_directory_structure()
    if dir_errors:
        all_passed = False
        print("   ❌ FAILED")
        for error in dir_errors:
            print(f"      - {error}")
    else:
        print("   ✅ PASSED")

    # Check Python syntax
    print("\n2. Checking Python syntax...")
    syntax_errors = check_python_syntax()
    if syntax_errors:
        all_passed = False
        print("   ❌ FAILED")
        for error in syntax_errors:
            print(f"      - {error}")
    else:
        print("   ✅ PASSED")

    # Check imports
    print("\n3. Checking module imports...")
    import_results = check_imports()
    import_failures = [r for r in import_results if not r[1]]

    for module, success, message in import_results:
        status = "✅" if success else "⚠️"
        print(f"   {status} {module}: {message}")

    # Only fail on unexpected import errors, not RDKit/numpy missing
    critical_failures = [
        r
        for r in import_failures
        if "rdkit" not in r[2].lower()
        and "numpy" not in r[2].lower()
        and "expected" not in r[2].lower()
    ]
    if critical_failures:
        all_passed = False
        print("   ❌ Critical import failures found")
    elif import_failures:
        print("   ⚠️  Some imports failed (missing dependencies - expected)")

    # Check test structure
    print("\n4. Checking test structure...")
    test_errors = check_test_structure()
    if test_errors:
        all_passed = False
        print("   ❌ FAILED")
        for error in test_errors:
            print(f"      - {error}")
    else:
        print("   ✅ PASSED")

    # Check required files
    print("\n5. Checking required files...")
    file_errors = check_required_files()
    if file_errors:
        all_passed = False
        print("   ❌ FAILED")
        for error in file_errors:
            print(f"      - {error}")
    else:
        print("   ✅ PASSED")

    # Summary
    print("\n" + "=" * 70)
    if all_passed:
        print("✅ ALL VALIDATIONS PASSED")
        print("\nProject structure is valid!")
        print("\nNext steps:")
        print("  1. Install dependencies: poetry install")
        print("  2. Run tests: make test")
        print("  3. Or use Docker: make docker-test")
        return 0
    else:
        print("❌ SOME VALIDATIONS FAILED")
        print("\nPlease fix the errors above before running tests.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
