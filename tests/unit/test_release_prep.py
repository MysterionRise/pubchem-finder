"""scripts/release_prep.py: roll the changelog, set every version, check a release."""

import importlib.util
import shutil
import sys
from pathlib import Path

import pytest

from pubchem_finder import __version__

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "release_prep.py"
FILES = ("CHANGELOG.md", "pyproject.toml", "server.json", "CITATION.cff")
URL = "https://github.com/MysterionRise/pubchem-finder"


def load_script():
    spec = importlib.util.spec_from_file_location("release_prep", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


rp = load_script()

CHANGELOG = f"""# Changelog

## [Unreleased]

### Fixed

- A fix.

## [1.2.0] - 2026-01-01

### Added

- A feature.

[Unreleased]: {URL}/compare/v1.2.0...HEAD
[1.2.0]: {URL}/releases/tag/v1.2.0
"""


def test_the_repository_is_ready_for_its_own_version():
    """The release workflow runs `check`; the tests run it on every PR."""
    assert rp.check(ROOT, __version__) == []


def test_roll_changelog_moves_entries_and_links():
    rolled = rp.roll_changelog(CHANGELOG, "1.2.1", "2026-02-03")
    assert "## [Unreleased]\n\n## [1.2.1] - 2026-02-03\n\n### Fixed\n" in rolled
    assert rp.changelog_section(rolled, "Unreleased") == ""
    assert rp.changelog_section(rolled, "1.2.1") == "### Fixed\n\n- A fix.\n"
    assert f"[Unreleased]: {URL}/compare/v1.2.1...HEAD\n" in rolled
    assert f"[1.2.1]: {URL}/releases/tag/v1.2.1\n[1.2.0]:" in rolled


def test_roll_changelog_refuses_empty_or_repeated_releases():
    empty = CHANGELOG.replace("### Fixed\n\n- A fix.\n\n", "")
    with pytest.raises(rp.ReleaseError, match="empty"):
        rp.roll_changelog(empty, "1.2.1", "2026-02-03")
    with pytest.raises(rp.ReleaseError, match="already has"):
        rp.roll_changelog(CHANGELOG, "1.2.0", "2026-02-03")


def test_changelog_section_matches_the_release_notes_rule():
    """Up to the next '## ' heading, like the release workflow's old awk rule."""
    assert rp.changelog_section(CHANGELOG, "Unreleased") == "### Fixed\n\n- A fix.\n"
    assert rp.changelog_section(CHANGELOG, "1.2.0").startswith(
        "### Added\n\n- A feature.\n"
    )
    assert rp.changelog_section(CHANGELOG, "9.9.9") == ""


def test_pyproject_version_only_in_the_project_table():
    text = '[project]\nname = "x"\nversion = "1.2.0"\n\n[tool.other]\nversion = "0"\n'
    updated = rp.set_pyproject_version(text, "1.2.1")
    assert 'version = "1.2.1"' in updated
    assert 'version = "0"' in updated
    assert rp.pyproject_version(updated) == "1.2.1"


def test_server_json_round_trips_unchanged():
    text = (ROOT / "server.json").read_text(encoding="utf-8")
    assert rp.set_server_json_version(text, __version__) == text
    moved = rp.set_server_json_version(text, "9.9.9")
    assert moved.count('"version": "9.9.9"') == text.count(
        f'"version": "{__version__}"'
    )


def test_citation_version_and_date():
    text = 'title: "x"\nversion: 1.2.0\ndate-released: "2026-01-01"\n'
    updated = rp.set_citation(text, "1.2.1", "2026-02-03")
    assert updated == 'title: "x"\nversion: 1.2.1\ndate-released: "2026-02-03"\n'


@pytest.fixture
def repo(tmp_path):
    for name in FILES:
        shutil.copy(ROOT / name, tmp_path / name)
    changelog = tmp_path / "CHANGELOG.md"
    text = changelog.read_text(encoding="utf-8")
    # Give the copy something to release.
    changelog.write_text(
        text.replace(
            "## [Unreleased]\n", "## [Unreleased]\n\n### Fixed\n\n- Next fix.\n", 1
        ),
        encoding="utf-8",
    )
    return tmp_path


def _next_patch(version: str) -> str:
    major, minor, patch = (int(p) for p in version.split("."))
    return f"{major}.{minor}.{patch + 1}"


def test_prep_then_check_passes(repo):
    version = _next_patch(__version__)
    calls = []
    rp.prep(repo, version, "2030-01-02", run=lambda cmd, **kw: calls.append((cmd, kw)))
    assert calls == [(["uv", "lock"], {"cwd": repo, "check": True})]
    notes = repo / "notes.md"
    assert rp.check(repo, version, notes) == []
    released = rp.changelog_section(
        (repo / "CHANGELOG.md").read_text(encoding="utf-8"), version
    )
    assert notes.read_text(encoding="utf-8") == released
    assert released.startswith("### Fixed\n\n- Next fix.\n")


def test_prep_refuses_old_or_malformed_versions_and_writes_nothing(repo):
    before = {name: (repo / name).read_text(encoding="utf-8") for name in FILES}
    with pytest.raises(rp.ReleaseError, match="not newer"):
        rp.prep(repo, __version__, "2030-01-02", lock=False)
    with pytest.raises(rp.ReleaseError, match="not a release version"):
        rp.prep(repo, "3.0", "2030-01-02", lock=False)
    with pytest.raises(rp.ReleaseError, match="not a date"):
        rp.prep(repo, _next_patch(__version__), "tomorrow", lock=False)
    assert before == {name: (repo / name).read_text(encoding="utf-8") for name in FILES}


def test_check_lists_every_problem(repo):
    problems = rp.check(repo, _next_patch(__version__))
    joined = "\n".join(problems)
    for name in ("pyproject.toml", "server.json", "CITATION.cff", "CHANGELOG.md"):
        assert name in joined
