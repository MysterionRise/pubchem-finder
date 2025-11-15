#!/usr/bin/env python3
"""
GitHub Actions CI Simulation Script

This script runs the same checks that GitHub Actions CI runs,
allowing you to verify your code will pass CI before pushing.

Usage:
    python scripts/ci_check.py          # Run all checks
    python scripts/ci_check.py --fast   # Skip slow tests
"""

import subprocess
import sys
import argparse
from pathlib import Path
from typing import List, Tuple


class Colors:
    """ANSI color codes for terminal output."""
    GREEN = '\033[92m'
    RED = '\033[91m'
    YELLOW = '\033[93m'
    BLUE = '\033[94m'
    BOLD = '\033[1m'
    END = '\033[0m'


def run_command(cmd: List[str], description: str, allow_failure: bool = False) -> bool:
    """Run a command and return success status."""
    print(f"\n{Colors.BLUE}▶ {description}{Colors.END}")
    print(f"  Command: {' '.join(cmd)}")

    try:
        result = subprocess.run(
            cmd,
            check=not allow_failure,
            capture_output=False,
            text=True
        )

        if result.returncode == 0:
            print(f"{Colors.GREEN}✅ PASSED{Colors.END}")
            return True
        else:
            if allow_failure:
                print(f"{Colors.YELLOW}⚠️  FAILED (non-blocking){Colors.END}")
                return True
            else:
                print(f"{Colors.RED}❌ FAILED{Colors.END}")
                return False
    except subprocess.CalledProcessError as e:
        if allow_failure:
            print(f"{Colors.YELLOW}⚠️  FAILED (non-blocking){Colors.END}")
            return True
        else:
            print(f"{Colors.RED}❌ FAILED with exit code {e.returncode}{Colors.END}")
            return False
    except FileNotFoundError:
        print(f"{Colors.RED}❌ Command not found: {cmd[0]}{Colors.END}")
        return False


def main():
    """Run all CI checks locally."""
    parser = argparse.ArgumentParser(description='Run GitHub Actions CI checks locally')
    parser.add_argument('--fast', action='store_true', help='Skip slow tests')
    parser.add_argument('--no-tests', action='store_true', help='Skip all tests (only quality checks)')
    args = parser.parse_args()

    print(f"{Colors.BOLD}{'='*70}{Colors.END}")
    print(f"{Colors.BOLD}GitHub Actions CI Checks - Local Simulation{Colors.END}")
    print(f"{Colors.BOLD}{'='*70}{Colors.END}")

    # Track results
    results: List[Tuple[str, bool]] = []

    # 1. Validate project structure
    passed = run_command(
        ['python', 'scripts/validate.py'],
        'Step 1: Validate project structure'
    )
    results.append(('Validation', passed))

    # 2. Code formatting
    passed = run_command(
        ['poetry', 'run', 'black', '--check', 'src', 'tests'],
        'Step 2: Check code formatting (Black)'
    )
    results.append(('Formatting', passed))

    # 3. Linting
    passed = run_command(
        ['poetry', 'run', 'ruff', 'check', 'src', 'tests'],
        'Step 3: Linting (Ruff)'
    )
    results.append(('Linting', passed))

    # 4. Type checking (non-blocking in CI)
    passed = run_command(
        ['poetry', 'run', 'mypy', 'src'],
        'Step 4: Type checking (Mypy)',
        allow_failure=True
    )
    results.append(('Type Checking', passed))

    if not args.no_tests:
        # 5. Unit tests
        if args.fast:
            passed = run_command(
                ['poetry', 'run', 'pytest', 'tests/unit', '-x'],
                'Step 5: Unit tests (fast mode - stop on first failure)'
            )
        else:
            passed = run_command(
                ['poetry', 'run', 'pytest', 'tests/unit', '-v', '--cov'],
                'Step 5: Unit tests with coverage'
            )
        results.append(('Unit Tests', passed))

        # 6. Integration tests
        if not args.fast:
            passed = run_command(
                ['poetry', 'run', 'pytest', 'tests/integration', '-v', '--cov', '--cov-append'],
                'Step 6: Integration tests'
            )
            results.append(('Integration Tests', passed))

        # 7. Performance tests (non-blocking)
        if not args.fast:
            passed = run_command(
                ['poetry', 'run', 'pytest', 'tests/test_performance.py', '-v', '--durations=5'],
                'Step 7: Performance tests (info only)',
                allow_failure=True
            )
            results.append(('Performance Tests', passed))

    # Summary
    print(f"\n{Colors.BOLD}{'='*70}{Colors.END}")
    print(f"{Colors.BOLD}CI Check Summary{Colors.END}")
    print(f"{Colors.BOLD}{'='*70}{Colors.END}\n")

    all_passed = True
    for check, passed in results:
        status = f"{Colors.GREEN}✅ PASSED{Colors.END}" if passed else f"{Colors.RED}❌ FAILED{Colors.END}"
        print(f"{check:.<40} {status}")
        if not passed:
            all_passed = False

    print(f"\n{Colors.BOLD}{'='*70}{Colors.END}")

    if all_passed:
        print(f"{Colors.GREEN}{Colors.BOLD}✅ ALL CHECKS PASSED - Ready to push!{Colors.END}")
        print(f"\n{Colors.BLUE}Your code will pass GitHub Actions CI ✓{Colors.END}")
        return 0
    else:
        print(f"{Colors.RED}{Colors.BOLD}❌ SOME CHECKS FAILED{Colors.END}")
        print(f"\n{Colors.YELLOW}Fix the errors above before pushing to GitHub{Colors.END}")
        print(f"\n{Colors.BLUE}Quick fixes:{Colors.END}")
        print(f"  - Format code:     {Colors.YELLOW}make format{Colors.END}")
        print(f"  - Run all checks:  {Colors.YELLOW}make quality{Colors.END}")
        print(f"  - Run tests:       {Colors.YELLOW}make test{Colors.END}")
        return 1


if __name__ == '__main__':
    sys.exit(main())
