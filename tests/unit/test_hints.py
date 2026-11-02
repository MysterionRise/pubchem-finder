"""Tests for the commands and costs the index messages give."""

from pathlib import Path

import pytest

from pubchem_finder import hints
from pubchem_finder.config import default_home


def test_release_commands_run_this_version_through_uvx(monkeypatch):
    monkeypatch.setattr(hints, "__version__", "2.2.0")
    assert (
        hints.command(hints.BUILD_SAMPLE)
        == "uvx --python 3.14 pubchem-finder==2.2.0 index build --limit 1 "
        "--with-substructure"
    )


@pytest.mark.parametrize("version", ["2.2.0.dev3+g1a2b3c", "0.0.0"])
def test_development_builds_run_the_installed_command(monkeypatch, version):
    monkeypatch.setattr(hints, "__version__", version)
    assert hints.command("index status") == "pubchem-finder index status"


def test_commands_name_the_data_directory_unless_it_is_the_default(monkeypatch):
    monkeypatch.setattr(hints, "__version__", "dev")
    assert (
        hints.command("index status", default_home()) == "pubchem-finder index status"
    )
    custom = Path("/data/pubchem")
    assert hints.command("index status", custom).endswith(f"--home {custom}")
    spaced = Path("/Users/me/My Data")
    assert hints.command("index status", spaced).endswith(f'--home "{spaced}"')


def test_build_hint_gives_the_command_and_its_cost():
    text = hints.build_index()
    assert hints.command(hints.BUILD_SAMPLE) in text
    assert "323 MB" in text
    assert "`--all`" in text


def test_launch_python_matches_the_published_commands():
    from tests.unit.test_metadata import LAUNCH_PYTHON

    assert hints.LAUNCH_PYTHON == LAUNCH_PYTHON
