"""Preview what ``--sample_metadata`` / ``--sample_pattern`` will resolve to.

A mirror of the matching rules in ``src/loading/sample_annotations.R``, kept
here so the form can populate its column picker and report the match rate the
moment you pick a file — a CSV header is one line and a regex is a regex, so
making you wait for a multi-minute workspace load to find out you typo'd a
column name would be absurd.

This is a *preview*, not the source of truth: the run itself always resolves
annotations in R. If the two ever disagree, R wins and this is the bug.
"""

from __future__ import annotations

import csv
import re
from pathlib import Path

#: Must match .ANNOT_UNASSIGNED in src/loading/sample_annotations.R.
UNASSIGNED = "Unassigned"

#: Column names treated as the sample identifier, in the R side's order.
ID_COLUMNS = ("sample_id", "sample", "sampleID", "sample_name", "filename",
              "file")

_PLACEHOLDERS = {"?", "na", "n/a", "null", "-", ""}


def annotation_key(value: str) -> str:
    """The join key: base name, no .fcs extension, no cytoset suffix, lowered.

    Mirrors ``annotation_key()`` in the R module.
    """
    value = re.sub(r".*[\\/]", "", str(value))
    value = re.sub(r"\.fcs_\d+$", "", value, flags=re.IGNORECASE)
    value = re.sub(r"\.fcs$", "", value, flags=re.IGNORECASE)
    return value.strip().lower()


def _read_header_and_keys(path: str) -> tuple[list[str], str, dict[str, dict]]:
    """Columns, the id column, and {key: row} for a metadata table."""
    text = Path(path).read_text(encoding="utf-8-sig", errors="replace")
    first = text.split("\n", 1)[0]
    delimiter = ";" if first.count(";") > first.count(",") else ","
    rows = list(csv.DictReader(text.splitlines(), delimiter=delimiter))
    if not rows:
        raise ValueError("the metadata file has a header but no rows")

    columns = [c for c in (rows[0].keys()) if c]
    id_column = next((c for c in ID_COLUMNS if c in columns), columns[0])
    keyed = {annotation_key(r.get(id_column, "")): r for r in reversed(rows)}
    return columns, id_column, keyed


def _clean(value) -> str:
    value = (value or "").strip()
    return UNASSIGNED if value.lower() in _PLACEHOLDERS else value


def describe(metadata: str, pattern: str, samples: list[str]) -> tuple[list, dict]:
    """``(columns, preview)`` for the current metadata/pattern settings.

    ``samples`` are the loaded sample display labels, which may be empty before
    a workspace has been read — the columns are still reported, only the match
    rate needs them. A CSV wins over a pattern, matching the R side.
    """
    if metadata:
        if not Path(metadata).is_file():
            return [], {"error": f"Metadata file does not exist: {metadata}"}
        try:
            columns, id_column, keyed = _read_header_and_keys(metadata)
        except (OSError, ValueError, csv.Error) as exc:
            return [], {"error": f"Could not read {Path(metadata).name}: {exc}"}

        value_cols = [c for c in columns if c != id_column]
        if not value_cols:
            return [], {"error": (
                f"{Path(metadata).name} has only the id column "
                f"('{id_column}') — there is nothing to annotate by.")}

        unmatched = [s for s in samples if annotation_key(s) not in keyed]
        levels = {
            c: len({_clean(r.get(c)) for r in keyed.values()})
            for c in value_cols
        }
        return value_cols, {
            "source": f"From {Path(metadata).name}",
            "total": len(samples),
            "matched": len(samples) - len(unmatched),
            "unmatched": unmatched,
            "levels": levels,
        }

    if pattern:
        # R (PCRE) writes named groups as (?<name>…); Python needs (?P<name>…),
        # and rejects the R form outright rather than just failing to name the
        # group. Both are accepted so the pattern you paste into the form is
        # the same one you can paste into the CLI.
        converted = re.sub(r"\(\?<(?=[A-Za-z_])", "(?P<", pattern)
        rx = None
        errors = []
        for candidate in (pattern, converted):
            try:
                compiled = re.compile(candidate)
            except re.error as exc:
                errors.append(str(exc))
                continue
            rx = compiled
            if compiled.groupindex:
                break
        if rx is None:
            return [], {"error": f"Not a valid regular expression: {errors[0]}"}
        if not rx.groupindex:
            return [], {"error": (
                "The pattern has no named groups, so there is nothing to name "
                "the columns after. Use (?<patient>…) style groups.")}

        columns = sorted(rx.groupindex, key=lambda k: rx.groupindex[k])
        hits = [rx.search(s) for s in samples]
        unmatched = [s for s, h in zip(samples, hits) if h is None]
        levels = {
            c: len({m.group(c) for m in hits if m and m.group(c) is not None})
            for c in columns
        }
        return columns, {
            "source": "From the sample name",
            "total": len(samples),
            "matched": len(samples) - len(unmatched),
            "unmatched": unmatched,
            "levels": levels,
        }

    return [], {}
