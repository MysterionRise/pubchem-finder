"""Data-directory resolution stays quiet; launchers own its diagnostic."""

import warnings
from pathlib import Path

import pytest

from pubchem_finder.config import default_home, relative_home_warning, resolve_home


@pytest.mark.parametrize("source", ["explicit", "environment"])
@pytest.mark.parametrize("warning_filter", ["default", "error", "ignore"])
def test_relative_home_resolves_quietly(
    source, warning_filter, tmp_path, monkeypatch, capsys
):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("PUBCHEM_FINDER_HOME", "data")
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter(warning_filter)
        for _ in range(2):
            assert (
                resolve_home("data" if source == "explicit" else None)
                == tmp_path / "data"
            )
    assert not caught
    assert capsys.readouterr() == ("", "")
    assert not (tmp_path / "data").exists()


def test_relative_home_diagnostic_has_no_output(tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    message = relative_home_warning("data")
    assert "relative data directory 'data'" in message
    assert str(tmp_path / "data") in message
    assert "working directory" in message
    assert "absolute path" in message
    assert capsys.readouterr() == ("", "")
    assert not (tmp_path / "data").exists()


@pytest.mark.parametrize("home", [None, "absolute", "tilde", ""])
def test_unambiguous_home_has_no_diagnostic(home, tmp_path, monkeypatch):
    monkeypatch.delenv("PUBCHEM_FINDER_HOME", raising=False)
    value = {"absolute": str(tmp_path), "tilde": "~/pubchem-test"}.get(home, home)
    expected = (
        default_home()
        if not value or not value.strip()
        else Path(value).expanduser().resolve()
    )
    assert resolve_home(value) == expected
    assert relative_home_warning(value) is None


def test_absolute_override_ignores_relative_environment(tmp_path, monkeypatch, capsys):
    monkeypatch.setenv("PUBCHEM_FINDER_HOME", "relative-env")
    assert resolve_home(tmp_path) == tmp_path
    assert relative_home_warning(tmp_path) is None
    assert capsys.readouterr() == ("", "")


def test_whitespace_environment_home_still_names_its_relative_directory(
    tmp_path, monkeypatch
):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("PUBCHEM_FINDER_HOME", "   ")
    expected = (tmp_path / "   ").resolve()
    assert resolve_home() == expected
    assert resolve_home() != default_home()
    assert str(expected) in relative_home_warning("   ")
