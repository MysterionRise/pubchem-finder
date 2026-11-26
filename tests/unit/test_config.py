"""Configured data directories must not silently depend on a client's CWD."""

import warnings
from pathlib import Path

import pytest

from pubchem_finder.config import default_home, resolve_home


@pytest.mark.parametrize("source", ["explicit", "environment"])
def test_relative_home_warns_and_keeps_cwd_resolution(source, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("PUBCHEM_FINDER_HOME", "data")
    with pytest.warns(UserWarning, match="relative data directory") as caught:
        resolved = resolve_home("data" if source == "explicit" else None)
    assert resolved == tmp_path / "data"
    message = str(caught[0].message)
    assert str(resolved) in message
    assert "working directory" in message
    assert "absolute path" in message
    assert not resolved.exists()


def test_relative_home_warning_is_not_repeated(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("default")
        for _ in range(2):
            resolve_home("data")
    assert len(caught) == 1


@pytest.mark.parametrize("home", [None, "absolute", "tilde"])
def test_unambiguous_home_does_not_warn(home, tmp_path, monkeypatch):
    monkeypatch.delenv("PUBCHEM_FINDER_HOME", raising=False)
    value = {None: None, "absolute": str(tmp_path), "tilde": "~/pubchem-test"}[home]
    expected = default_home() if value is None else Path(value).expanduser().resolve()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        assert resolve_home(value) == expected
    assert not caught


def test_absolute_override_does_not_warn_about_ignored_environment(
    tmp_path, monkeypatch
):
    monkeypatch.setenv("PUBCHEM_FINDER_HOME", "relative-env")
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        assert resolve_home(tmp_path) == tmp_path
    assert not caught
