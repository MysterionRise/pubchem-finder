"""Prepare a release, or check that one is ready.

    make release-prep VERSION=2.1.2
    uv run --no-project python scripts/release_prep.py check 2.1.2 --notes release-notes.md

``prep`` moves CHANGELOG.md's Unreleased entries under ``## [x.y.z] - <date>`` and
updates the compare links, sets the version in pyproject.toml, server.json (every
version field) and CITATION.cff (with ``date-released``), then runs ``uv lock``.
``check`` verifies all of that and writes the release notes; the release workflow runs
it, and the tests run it on the current version. The Claude plugin moves separately,
once the release is on PyPI: ``make plugin-lock VERSION=x.y.z``.

Standard library only, so it runs with ``uv run --no-project``.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import re
import subprocess
import sys
from collections.abc import Callable
from pathlib import Path
from typing import Any

REPO_URL = "https://github.com/MysterionRise/pubchem-finder"
ROOT = Path(__file__).resolve().parents[1]
_VERSION = re.compile(r"\d+\.\d+\.\d+")
_DATE = re.compile(r"\d{4}-\d{2}-\d{2}")


class ReleaseError(Exception):
    """The files aren't in a state a release can be made from."""


def parse_version(version: str) -> tuple[int, int, int]:
    if not _VERSION.fullmatch(version):
        raise ReleaseError(f"not a release version: {version!r} (want x.y.z)")
    major, minor, patch = (int(part) for part in version.split("."))
    return major, minor, patch


def changelog_section(text: str, version: str) -> str:
    """The lines under ``## [version]`` up to the next ``## `` heading, trimmed."""
    lines: list[str] = []
    inside = False
    for line in text.splitlines():
        if line.startswith(f"## [{version}]"):
            inside = True
            continue
        if inside and line.startswith("## "):
            break
        if inside:
            lines.append(line)
    body = "\n".join(lines).strip("\n")
    return body + "\n" if body else ""


def changelog_date(text: str, version: str) -> str | None:
    match = re.search(
        rf"^## \[{re.escape(version)}\] - (\S+)$", text, flags=re.MULTILINE
    )
    return match.group(1) if match else None


def roll_changelog(text: str, version: str, date: str) -> str:
    """Move the Unreleased entries into a new ``## [version] - date`` section."""
    heading = "## [Unreleased]\n"
    if text.count(heading) != 1:
        raise ReleaseError("CHANGELOG.md needs exactly one '## [Unreleased]' heading")
    if f"## [{version}]" in text:
        raise ReleaseError(f"CHANGELOG.md already has a section for {version}")
    if not changelog_section(text, "Unreleased").strip():
        raise ReleaseError("CHANGELOG.md's Unreleased section is empty")
    link = re.search(
        r"^\[Unreleased\]: \S+/compare/v\d+\.\d+\.\d+\.\.\.HEAD$",
        text,
        flags=re.MULTILINE,
    )
    if link is None:
        raise ReleaseError(
            "CHANGELOG.md has no '[Unreleased]: .../compare/vX...HEAD' link"
        )
    text = text.replace(heading, f"{heading}\n## [{version}] - {date}\n", 1)
    return text.replace(
        link.group(0),
        f"[Unreleased]: {REPO_URL}/compare/v{version}...HEAD\n"
        f"[{version}]: {REPO_URL}/releases/tag/v{version}",
        1,
    )


def _project_table(text: str) -> tuple[int, int]:
    start = re.search(r"^\[project\]\s*$", text, flags=re.MULTILINE)
    if start is None:
        raise ReleaseError("pyproject.toml has no [project] table")
    end = re.search(r"^\[", text[start.end() :], flags=re.MULTILINE)
    return start.end(), start.end() + end.start() if end else len(text)


def pyproject_version(text: str) -> str:
    lo, hi = _project_table(text)
    found = re.findall(r'^version = "([^"]*)"$', text[lo:hi], flags=re.MULTILINE)
    if len(found) != 1:
        raise ReleaseError("pyproject.toml's [project] table needs exactly one version")
    return str(found[0])


def set_pyproject_version(text: str, version: str) -> str:
    lo, hi = _project_table(text)
    table, count = re.subn(
        r'^version = "[^"]*"$',
        f'version = "{version}"',
        text[lo:hi],
        flags=re.MULTILINE,
    )
    if count != 1:
        raise ReleaseError("pyproject.toml's [project] table needs exactly one version")
    return text[:lo] + table + text[hi:]


def _server_versions(data: dict[str, Any]) -> dict[str, str]:
    return {
        "server.json": data["version"],
        **{
            f"server.json packages[{i}]": package["version"]
            for i, package in enumerate(data["packages"])
        },
    }


def set_server_json_version(text: str, version: str) -> str:
    data = json.loads(text)
    data["version"] = version
    for package in data["packages"]:
        package["version"] = version
    return json.dumps(data, indent=2, ensure_ascii=False) + "\n"


def _citation_field(text: str, field: str) -> str | None:
    match = re.search(rf'^{field}: "?([^"\n]*)"?$', text, flags=re.MULTILINE)
    return match.group(1) if match else None


def set_citation(text: str, version: str, date: str) -> str:
    text, n_version = re.subn(
        r"^version: .*$", f"version: {version}", text, flags=re.MULTILINE
    )
    text, n_date = re.subn(
        r"^date-released: .*$", f'date-released: "{date}"', text, flags=re.MULTILINE
    )
    if (n_version, n_date) != (1, 1):
        raise ReleaseError("CITATION.cff needs one version and one date-released line")
    return text


def prep(
    root: Path,
    version: str,
    date: str,
    *,
    lock: bool = True,
    run: Callable[..., Any] = subprocess.run,
) -> None:
    """Roll the changelog and set every version; nothing is written on an error."""
    if not _DATE.fullmatch(date):
        raise ReleaseError(f"not a date: {date!r} (want YYYY-MM-DD)")
    files = {
        name: (root / name).read_text(encoding="utf-8")
        for name in ("CHANGELOG.md", "pyproject.toml", "server.json", "CITATION.cff")
    }
    current = pyproject_version(files["pyproject.toml"])
    if parse_version(version) <= parse_version(current):
        raise ReleaseError(f"{version} is not newer than the current version {current}")
    updated = {
        "CHANGELOG.md": roll_changelog(files["CHANGELOG.md"], version, date),
        "pyproject.toml": set_pyproject_version(files["pyproject.toml"], version),
        "server.json": set_server_json_version(files["server.json"], version),
        "CITATION.cff": set_citation(files["CITATION.cff"], version, date),
    }
    for name, text in updated.items():
        (root / name).write_text(text, encoding="utf-8")
    if lock:
        run(["uv", "lock"], cwd=root, check=True)


def check(root: Path, version: str, notes: Path | None = None) -> list[str]:
    """Problems that stop ``version`` from being released (empty when ready)."""
    parse_version(version)
    problems: list[str] = []
    found = {
        "pyproject.toml": pyproject_version(
            (root / "pyproject.toml").read_text(encoding="utf-8")
        ),
        **_server_versions(
            json.loads((root / "server.json").read_text(encoding="utf-8"))
        ),
    }
    citation = (root / "CITATION.cff").read_text(encoding="utf-8")
    found["CITATION.cff"] = _citation_field(citation, "version") or ""
    problems += [
        f"{name} says {v}, not {version}" for name, v in found.items() if v != version
    ]

    changelog = (root / "CHANGELOG.md").read_text(encoding="utf-8")
    section = changelog_section(changelog, version)
    date = changelog_date(changelog, version)
    if not section.strip():
        problems.append(f"CHANGELOG.md has no entries under '## [{version}]'")
    if date is None or not _DATE.fullmatch(date):
        problems.append(f"CHANGELOG.md needs a '## [{version}] - YYYY-MM-DD' heading")
    released = _citation_field(citation, "date-released")
    if released != date:
        problems.append(
            f"CITATION.cff date-released is {released}, the changelog says {date}"
        )
    if f"[{version}]: {REPO_URL}/releases/tag/v{version}" not in changelog:
        problems.append(f"CHANGELOG.md has no [{version}] release link")
    if f"[Unreleased]: {REPO_URL}/compare/v{version}...HEAD" not in changelog:
        problems.append(
            f"CHANGELOG.md's [Unreleased] link doesn't compare from v{version}"
        )
    if notes is not None and section.strip():
        notes.write_text(section, encoding="utf-8")
    return problems


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    commands = parser.add_subparsers(dest="command", required=True)
    prep_cmd = commands.add_parser("prep", help="roll the changelog, set every version")
    prep_cmd.add_argument("version")
    prep_cmd.add_argument("--date", default=dt.date.today().isoformat())
    prep_cmd.add_argument("--no-lock", action="store_true", help="skip `uv lock`")
    check_cmd = commands.add_parser("check", help="verify a release is ready")
    check_cmd.add_argument("version")
    check_cmd.add_argument("--notes", type=Path, help="write the release notes here")
    args = parser.parse_args(argv)
    try:
        if args.command == "prep":
            prep(ROOT, args.version, args.date, lock=not args.no_lock)
            print(
                f"ready for {args.version}: review the diff, then open the release PR"
            )
            return 0
        problems = check(ROOT, args.version, args.notes)
    except ReleaseError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    for problem in problems:
        print(f"::error::{problem}")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
