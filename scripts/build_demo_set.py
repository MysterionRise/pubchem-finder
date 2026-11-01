"""Regenerate the demo set: FDA-approved drugs from PubChem, with their structures.

    uv run python scripts/build_demo_set.py

Takes the compounds PubChem links to Drugs@FDA under its "FDA Approved Drugs"
heading and replaces each with its PubChem parent compound (salts and solvates
removed: metformin, not metformin hydrochloride), since that is what people search
for. Single-component compounds PubChem gives no parent (zwitterions such as
carnitine, elements) stay as listed; mixtures without a parent are left out.
Writes PubChem's SMILES, InChIKey, IUPAC name, formula and molecular weight to
``src/pubchem_finder/index/demo_data/`` as a TSV (one row per CID, sorted) and a
JSON file saying where it came from. Drugs@FDA is a U.S. government work in the
public domain. Run it before a minor release and review the diff: the demo index
is rebuilt on each user's machine when the TSV changes.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import io
import json
import sys
from pathlib import Path
from typing import Any

from pubchem_finder.data.pug_rest import PubChemClient
from pubchem_finder.index.builder import METADATA_SCHEMA

OUT = Path(__file__).resolve().parents[1] / "src/pubchem_finder/index/demo_data"
NAME = "fda_approved_drugs"
HEADING = "FDA Approved Drugs"
SOURCE = "Drugs@FDA"
PROPERTIES = ("SMILES", "InChIKey", "IUPACName", "MolecularFormula", "MolecularWeight")
BATCH = 1000


def fda_cids(client: PubChemClient) -> tuple[list[int], dict[str, Any]]:
    """CIDs PubChem links to Drugs@FDA annotations, and the source's license."""
    url = f"{client.view_url}/annotations/heading/JSON"
    cids: set[int] = set()
    license_info: dict[str, Any] = {}
    page, pages = 1, 1
    while page <= pages:
        payload = client.request("GET", url, params={"heading": HEADING, "page": page})
        annotations = payload["Annotations"]
        pages = int(annotations.get("TotalPages", 1))
        for item in annotations["Annotation"]:
            if item.get("SourceName") != SOURCE:
                continue
            cids.update(int(c) for c in item.get("LinkedRecords", {}).get("CID", []))
            license_info = license_info or {
                "license": item.get("LicenseNote"),
                "license_url": item.get("LicenseURL"),
                "source_url": item.get("URL"),
            }
        page += 1
    return sorted(cids), license_info


def parent_cids(client: PubChemClient, cids: list[int]) -> set[int]:
    """PubChem's parent compounds of ``cids`` (a compound that is its own parent
    maps to itself; one without a parent maps to nothing)."""
    url = f"{client.base_url}/compound/cid/cids/JSON"
    parents: set[int] = set()
    for batch in _batches(cids):
        payload = client.request(
            "POST", url, params={"cids_type": "parent"}, data={"cid": batch}
        )
        parents.update(int(c) for c in payload["IdentifierList"]["CID"])
    return parents


def properties(client: PubChemClient, cids: list[int]) -> dict[int, dict[str, Any]]:
    """PubChem's properties for ``cids``, by CID."""
    url = f"{client.base_url}/compound/cid/property/{','.join(PROPERTIES)}/JSON"
    rows: dict[int, dict[str, Any]] = {}
    for batch in _batches(cids):
        payload = client.request("POST", url, data={"cid": batch})
        rows.update((int(r["CID"]), r) for r in payload["PropertyTable"]["Properties"])
    return rows


def _batches(cids: list[int]) -> list[str]:
    return [
        ",".join(str(c) for c in cids[start : start + BATCH])
        for start in range(0, len(cids), BATCH)
    ]


def demo_set(
    listed: list[int], parents: set[int], rows: dict[int, dict[str, Any]]
) -> list[dict[str, Any]]:
    """Parents, plus listed single-component compounds without one."""
    keep = set(parents)
    for cid in listed:
        smiles = rows.get(cid, {}).get("SMILES") or ""
        if cid not in parents and smiles and "." not in smiles:
            keep.add(cid)
    return [rows[cid] for cid in sorted(keep) if cid in rows]


def to_tsv(rows: list[dict[str, Any]]) -> str:
    """One row per CID, sorted, in the index's metadata columns."""
    buf = io.StringIO()
    writer = csv.writer(buf, delimiter="\t", lineterminator="\n")
    writer.writerow(METADATA_SCHEMA.names)
    seen: set[int] = set()
    for row in sorted(rows, key=lambda r: int(r["CID"])):
        cid = int(row["CID"])
        if cid in seen or not row.get("SMILES"):
            continue
        seen.add(cid)
        writer.writerow([cid, *(str(row.get(p) or "").strip() for p in PROPERTIES)])
    return buf.getvalue()


def build(client: PubChemClient, out: Path, today: str) -> int:
    """Write ``<NAME>.tsv`` and ``<NAME>.json`` to ``out``; returns the row count."""
    listed, license_info = fda_cids(client)
    parents = parent_cids(client, listed)
    rows = properties(client, sorted(set(listed) | parents))
    tsv = to_tsv(demo_set(listed, parents, rows))
    n_rows = tsv.count("\n") - 1
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{NAME}.tsv").write_bytes(tsv.encode("utf-8"))
    provenance = {
        "name": "FDA-approved drugs",
        "description": (
            f"PubChem parent compounds of the {len(listed)} compounds PubChem links "
            f"to {SOURCE} under its '{HEADING}' heading (single-component "
            "compounds without a parent kept as listed), with structures and "
            "properties computed by PubChem"
        ),
        "source": f"{SOURCE}, U.S. Food and Drug Administration, via PubChem",
        **license_info,
        "pubchem_heading": HEADING,
        "retrieved": today,
        "n_compounds": n_rows,
        "generator": "scripts/build_demo_set.py",
    }
    (out / f"{NAME}.json").write_text(
        json.dumps(provenance, indent=2) + "\n", encoding="utf-8"
    )
    return n_rows


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--out", type=Path, default=OUT, help="output directory")
    args = parser.parse_args(argv)
    with PubChemClient(timeout=120) as client:
        n_rows = build(client, args.out, dt.date.today().isoformat())
    print(f"{n_rows} compounds written to {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
