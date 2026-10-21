"""Pure functions that turn PubChem JSON into this package's models.

Kept free of HTTP so every parsing rule can be tested on recorded payloads.
"""

import re
from collections import Counter, defaultdict
from collections.abc import Iterable, Sequence
from typing import Any

from ..models import (
    Annotation,
    BioactivitySummary,
    PotentResult,
    SourcedText,
    TargetActivity,
)

CAS_RE = re.compile(r"^(\d{2,7})-(\d{2})-(\d)$")


class PayloadError(ValueError):
    """The payload does not have the expected structure."""


def _information(payload: dict[str, Any]) -> list[dict[str, Any]]:
    try:
        info = payload["InformationList"]["Information"]
    except (KeyError, TypeError) as exc:
        raise PayloadError(f"missing InformationList: {exc}") from exc
    if not isinstance(info, list):
        raise PayloadError("InformationList.Information is not a list")
    return info


def is_cas_number(text: str) -> bool:
    """True if ``text`` is a CAS Registry Number with a valid check digit."""
    match = CAS_RE.match(text.strip())
    if not match:
        return False
    digits = (match.group(1) + match.group(2))[::-1]
    total = sum((i + 1) * int(d) for i, d in enumerate(digits))
    return total % 10 == int(match.group(3))


def parse_synonyms(payload: dict[str, Any]) -> list[str]:
    """Synonyms in PubChem's order, case-insensitive duplicates removed."""
    seen: set[str] = set()
    result: list[str] = []
    for entry in _information(payload):
        for name in entry.get("Synonym", []):
            text = str(name).strip()
            if text and text.lower() not in seen:
                seen.add(text.lower())
                result.append(text)
    return result


def parse_xrefs(payload: dict[str, Any], types: Sequence[str]) -> dict[str, list[str]]:
    """Identifiers per xref type, deduplicated, order preserved."""
    out: dict[str, list[str]] = {t: [] for t in types}
    for entry in _information(payload):
        for t in types:
            seen = set(out[t])
            for value in entry.get(t, []):
                text = str(value)
                if text not in seen:
                    seen.add(text)
                    out[t].append(text)
    return out


def parse_descriptions(
    payload: dict[str, Any], max_items: int = 5, max_chars: int = 1500
) -> tuple[str | None, list[SourcedText]]:
    """Compound title and up to ``max_items`` sourced descriptions."""
    title: str | None = None
    texts: list[SourcedText] = []
    for entry in _information(payload):
        if title is None and entry.get("Title"):
            title = str(entry["Title"])
        text = entry.get("Description")
        if text and len(texts) < max_items:
            body = str(text).strip()
            if len(body) > max_chars:
                body = body[: max_chars - 1].rstrip() + "…"
            texts.append(
                SourcedText(
                    text=body,
                    source=entry.get("DescriptionSourceName"),
                    url=entry.get("DescriptionURL"),
                )
            )
    return title, texts


# --------------------------------------------------------------------------- #
# Bioassay summary
# --------------------------------------------------------------------------- #

_COLUMN_ALIASES: dict[str, tuple[str, ...]] = {
    "aid": ("AID",),
    "outcome": ("Activity Outcome", "Bioactivity Outcome"),
    "accession": ("Target Accession", "Protein Accession"),
    "gene_id": ("Target GeneID", "Gene ID"),
    "value": ("Activity Value [uM]", "Activity Value [µM]"),
    "activity_name": ("Activity Name",),
    "assay_name": ("Assay Name", "BioAssay Name"),
}


def _column_index(columns: Sequence[str]) -> dict[str, int]:
    lowered = {c.strip().lower(): i for i, c in enumerate(columns)}
    index: dict[str, int] = {}
    for key, names in _COLUMN_ALIASES.items():
        for name in names:
            if name.lower() in lowered:
                index[key] = lowered[name.lower()]
                break
    return index


def summarize_assays(
    payload: dict[str, Any], cid: int, top_n: int = 5
) -> BioactivitySummary:
    """Aggregate PubChem's ``assaysummary`` table for one CID.

    The raw table has one row per (assay, substance) and can run to thousands
    of rows for well-studied compounds; the summary keeps what an assistant
    needs: how often the compound was tested, outcome counts, the targets with
    the most active results, and the most potent active measurements.
    """
    try:
        table = payload["Table"]
        columns: list[str] = list(table["Columns"]["Column"])
        rows: list[dict[str, Any]] = list(table.get("Row", []))
    except (KeyError, TypeError) as exc:
        raise PayloadError(f"missing assay Table: {exc}") from exc
    col = _column_index(columns)
    if "aid" not in col or "outcome" not in col:
        raise PayloadError(f"assay table lacks AID/outcome columns: {columns}")

    def cell(cells: Sequence[Any], key: str) -> str | None:
        i = col.get(key)
        if i is None or i >= len(cells):
            return None
        value = str(cells[i]).strip()
        return value or None

    aids: set[int] = set()
    outcomes: Counter[str] = Counter()
    target_active: Counter[tuple[str | None, str | None]] = Counter()
    target_aids: dict[tuple[str | None, str | None], set[int]] = defaultdict(set)
    potent: list[PotentResult] = []
    for row in rows:
        cells = row.get("Cell", [])
        aid_text = cell(cells, "aid")
        if not aid_text or not aid_text.isdigit():
            continue
        aid = int(aid_text)
        aids.add(aid)
        outcome = (cell(cells, "outcome") or "Unspecified").capitalize()
        outcomes[outcome] += 1
        if outcome != "Active":
            continue
        target = (cell(cells, "accession"), cell(cells, "gene_id"))
        if any(target) and aid not in target_aids[target]:
            target_aids[target].add(aid)
            target_active[target] += 1
        value_text = cell(cells, "value")
        try:
            value = float(value_text) if value_text else None
        except ValueError:
            value = None
        if value is not None:
            potent.append(
                PotentResult(
                    aid=aid,
                    assay_name=cell(cells, "assay_name"),
                    activity_name=cell(cells, "activity_name"),
                    value_um=value,
                    target_accession=target[0],
                )
            )
    potent.sort(key=lambda p: p.value_um)
    return BioactivitySummary(
        cid=cid,
        tested=bool(aids),
        n_assays=len(aids),
        outcomes=dict(outcomes.most_common()),
        top_targets=[
            TargetActivity(target_accession=acc, target_gene_id=gene, active_assays=n)
            for (acc, gene), n in target_active.most_common(top_n)
        ],
        most_potent=potent[:top_n],
    )


# --------------------------------------------------------------------------- #
# PUG View
# --------------------------------------------------------------------------- #


def _value_text(value: dict[str, Any]) -> str:
    parts: list[str] = []
    for item in value.get("StringWithMarkup", []):
        text = str(item.get("String", "")).strip()
        if not text:
            # Pictograms carry their meaning in the markup, not the string.
            extras = [str(m["Extra"]) for m in item.get("Markup", []) if m.get("Extra")]
            text = ", ".join(extras)
        if text:
            parts.append(text)
    numbers = value.get("Number")
    if numbers:
        unit = str(value.get("Unit", "")).strip()
        joined = ", ".join(
            f"{n:g}" if isinstance(n, float) else str(n) for n in numbers
        )
        parts.append(f"{joined} {unit}".strip())
    return "; ".join(parts)


def parse_annotations(payload: dict[str, Any]) -> tuple[str | None, list[Annotation]]:
    """Flatten a PUG View record into annotations with their sources."""
    try:
        record = payload["Record"]
    except (KeyError, TypeError) as exc:
        raise PayloadError(f"missing Record: {exc}") from exc
    refs = {
        ref.get("ReferenceNumber"): ref
        for ref in record.get("Reference", [])
        if isinstance(ref, dict)
    }
    out: list[Annotation] = []

    def walk(sections: Iterable[dict[str, Any]], path: list[str]) -> None:
        for section in sections:
            here = [*path, str(section.get("TOCHeading", "")).strip()]
            for info in section.get("Information", []):
                text = _value_text(info.get("Value", {}))
                if not text:
                    continue
                ref = refs.get(info.get("ReferenceNumber"), {})
                out.append(
                    Annotation(
                        section=" > ".join(p for p in here if p),
                        name=info.get("Name"),
                        text=text,
                        source=ref.get("SourceName"),
                        url=ref.get("URL"),
                        license=ref.get("LicenseURL") or ref.get("LicenseNote"),
                    )
                )
            walk(section.get("Section", []), here)

    walk(record.get("Section", []), [])
    title = record.get("RecordTitle")
    return (str(title) if title else None), out


def parse_headings(payload: dict[str, Any]) -> list[str]:
    """All section headings in a PUG View index, in page order, deduplicated."""
    try:
        record = payload["Record"]
    except (KeyError, TypeError) as exc:
        raise PayloadError(f"missing Record: {exc}") from exc
    seen: set[str] = set()
    headings: list[str] = []

    def walk(sections: Iterable[dict[str, Any]]) -> None:
        for section in sections:
            heading = str(section.get("TOCHeading", "")).strip()
            if heading and heading not in seen:
                seen.add(heading)
                headings.append(heading)
            walk(section.get("Section", []))

    walk(record.get("Section", []))
    return headings


def trim_annotations(
    annotations: list[Annotation], max_chars: int
) -> tuple[list[Annotation], bool]:
    """Keep whole annotations, in order, until ``max_chars`` of text is used."""
    kept: list[Annotation] = []
    used = 0
    for annotation in annotations:
        size = len(annotation.text) + len(annotation.section) + 20
        if used + size > max_chars:
            if kept:
                return kept, True
            # A single oversized annotation: keep its beginning.
            room = max(200, max_chars - len(annotation.section) - 20)
            short = annotation.text[: room - 1].rstrip() + "…"
            return [annotation.model_copy(update={"text": short})], True
        kept.append(annotation)
        used += size
    return kept, False
