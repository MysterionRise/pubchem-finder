#!/bin/bash
# Quick CI status check script
# Shows what will run in GitHub Actions CI

set -e

echo "╔══════════════════════════════════════════════════════════════════╗"
echo "║           GitHub Actions CI - Configuration Check                ║"
echo "╚══════════════════════════════════════════════════════════════════╝"
echo ""

# Check if workflow file exists
if [ ! -f .github/workflows/ci.yml ]; then
    echo "❌ CI workflow file not found!"
    exit 1
fi

echo "✅ CI workflow file found: .github/workflows/ci.yml"
echo ""

# Show CI configuration
echo "📋 CI Configuration:"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

# Extract key info from CI file
echo "Triggers:"
grep -A 2 "^on:" .github/workflows/ci.yml | grep "branches:" -A 1 | grep -v "branches:" || true

echo ""
echo "Python Versions:"
grep "python-version:" .github/workflows/ci.yml | head -2

echo ""
echo "📊 Test Steps That Will Run:"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

# Extract all step names
grep "name:" .github/workflows/ci.yml | grep -v "^name:" | sed 's/.*name://' | sed 's/^[[:space:]]*//' | nl -w2 -s'. '

echo ""
echo "🎯 To simulate CI locally:"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  python scripts/ci_check.py           # Full simulation"
echo "  python scripts/ci_check.py --fast    # Quick check"
echo "  make quality && make test            # All checks"
echo ""

echo "✅ GitHub Actions CI is properly configured!"
echo ""
echo "📖 Documentation:"
echo "   - CI Details:      .github/workflows/README.md"
echo "   - Testing Guide:   TESTING.md"
echo "   - Verification:    CI_VERIFICATION.md"
