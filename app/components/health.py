"""Render the R load report: what came in, what did not, and what to use.

Presentation only. Every number, verdict and diagnosis shown here was computed
by ``src/visualization/umap_load_report.R``; this module decides what is worth
seeing first and nothing else. If you find yourself working something out in
here, it belongs in R — the app cannot see a GatingSet, so anything it derives
is a second opinion that will eventually disagree with the first.

The panel is built around one idea: **a glance, then a click**. The headline is
one line. Everything else — the workspace list, the marker panel, the sample
names behind each verdict — is behind its own disclosure, because a
thirty-three-workspace load has more detail than anyone wants at once and the
text log it replaces failed precisely by showing all of it.

Workspaces are switchable. Clicking one off excludes it from the run and from
the next load (see select_workspaces() in src/loading/flowjo_workspace.R), which
is the only way to act on "this workspace is junk" without editing a path.
"""

from __future__ import annotations

import dash_mantine_components as dmc
from dash import html

#: status -> (colour, one word). Deliberately one word: the row already carries
#: the counts, and a sentence per workspace is what made this unreadable.
_STATUS = {
    "loaded": ("green", "loaded"),
    "build_failed": ("red", "refused"),
    "unusable": ("red", "unreadable"),
    "no_files": ("orange", "no files"),
    "not_loaded": ("gray", "skipped"),
}

#: Sample names shown per verdict before the list is cut short.
_SAMPLE_LIMIT = 40


def _sample_list(title: str, names: list[str], color: str):
    """A collapsible list of sample names, or nothing when there are none."""
    if not names:
        return None
    shown = names[:_SAMPLE_LIMIT]
    extra = len(names) - len(shown)
    return dmc.Accordion(
        chevronPosition="left", variant="filled",
        children=[
            dmc.AccordionItem(
                value=title,
                children=[
                    dmc.AccordionControl(
                        dmc.Text(f"{title} ({len(names)})", size="xs", c=color)),
                    dmc.AccordionPanel(
                        dmc.Stack(
                            [dmc.Text(n, size="xs", ff="monospace")
                             for n in shown]
                            + ([dmc.Text(f"… {extra} more — see the Log tab",
                                         size="xs", c="dimmed", fs="italic")]
                               if extra else []),
                            gap=2)),
                ],
            )
        ],
    )


def _workspace_row(ws: dict, excluded: bool):
    """One workspace: a switch, its name, its counts, and detail on demand."""
    color, meaning = _STATUS.get(ws.get("status", ""), ("gray", ws.get("status", "")))
    routes = ws.get("routes") or {}
    loaded = ws.get("status") == "loaded"

    # Only the routes that contributed. A row of zeroes is noise.
    pills = [
        dmc.Badge(f"{routes[k]} {lbl}", color=c, variant="light", size="xs")
        for k, lbl, c in (("name", "by name", "green"),
                          ("recorded", "by recorded name", "teal"),
                          ("content", "by $FIL", "blue"))
        if routes.get(k)
    ]
    if ws.get("n_gates"):
        pills.append(dmc.Badge(f"{ws['n_gates']} gates", variant="light",
                               size="xs", color="gray"))
    if ws.get("markers"):
        pills.append(dmc.Badge(f"{len(ws['markers'])} markers", variant="light",
                               size="xs", color="gray"))
    # A workspace whose samples were not all acquired with the same channels.
    # It loaded — but as several gating sets, and that is a property of the data
    # worth seeing next to the counts rather than only in the log.
    if (ws.get("n_panel_parts") or 1) > 1:
        pills.append(dmc.Badge(f"{ws['n_panel_parts']} channel sets",
                               variant="light", size="xs", color="violet"))

    head = dmc.Group(
        [
            # A workspace that did not load cannot be included by ticking it, so
            # its switch is disabled rather than lying about being available.
            dmc.Switch(
                id={"comp": "wsp-toggle", "bid": ws["label"]},
                checked=loaded and not excluded,
                disabled=not loaded,
                size="xs",
            ),
            dmc.Text(ws["label"], fw=600, size="sm",
                     td="line-through" if (excluded or not loaded) else None,
                     c="dimmed" if (excluded or not loaded) else None),
            dmc.Badge(meaning, color=color, variant="light", size="xs"),
            dmc.Text(f"{ws.get('n_loaded') or 0}/{ws.get('n_listed') or 0}",
                     size="xs", c="dimmed"),
        ],
        gap="xs", wrap="nowrap",
    )

    detail = []
    if ws.get("reason"):
        detail.append(dmc.Text(ws["reason"], size="xs", c="red",
                               style={"whiteSpace": "pre-wrap"}))
    if ws.get("markers"):
        detail.append(dmc.Text("Panel: " + ", ".join(ws["markers"]), size="xs",
                               c="dimmed"))
    if ws.get("markers_not_shared"):
        detail.append(dmc.Text("Not shared: "
                               + ", ".join(ws["markers_not_shared"]),
                               size="xs", c="orange"))
    for title, key, colour in (
        ("Paired by content, not by name", "rescued", "blue"),
        ("Listed but did not load", "rejected", "red"),
        ("Loaded but never listed", "unlisted", "orange"),
    ):
        block = _sample_list(title, ws.get(key) or [], colour)
        if block is not None:
            detail.append(block)

    body = [head]
    if pills:
        body.append(dmc.Group(pills, gap=4, wrap="wrap"))
    body.extend(detail)
    return dmc.Paper(dmc.Stack(body, gap=6), p="xs", withBorder=True,
                     radius="sm",
                     style={"opacity": 0.55} if excluded else None)


def headline(inspect: dict | None, excluded: list[str] | None = None):
    """The one line that is always visible, above the collapsed sections."""
    report = (inspect or {}).get("report")
    if not report:
        return None
    totals = report.get("totals") or {}
    found = totals.get("workspaces_found") or 0
    loaded = totals.get("workspaces_loaded") or 0
    n_off = len([w for w in (excluded or [])])
    samples = totals.get("samples") or 0

    bits = [f"{loaded} of {found} workspaces", f"{samples} samples"]
    if n_off:
        bits.append(f"{n_off} switched off")
    colour = "green" if loaded == found and not n_off else (
        "red" if not loaded else "yellow")
    return dmc.Badge(" · ".join(bits), color=colour, variant="light", size="sm")


def marker_section(report: dict):
    """Common panel, what diverges, and the aliases R suggests for it."""
    markers = report.get("markers") or {}
    common = markers.get("common") or []
    not_common = markers.get("not_common") or []
    suggestions = markers.get("suggestions") or []
    conflicts = markers.get("conflicts") or []

    body = [
        dmc.Text(", ".join(common) if common
                 else "No marker is shared by every workspace — a joint UMAP "
                      "is not possible until this is resolved.",
                 size="xs", c="dimmed" if common else "red"),
    ]

    if not_common:
        body.append(dmc.Space(h=4))
        body.append(dmc.Stack(
            [dmc.Text(f"{m['marker']} — only in {len(m['present_in'])}: "
                      + ", ".join(m["present_in"][:4])
                      + ("…" if len(m["present_in"]) > 4 else ""),
                      size="xs", c="orange")
             for m in not_common],
            gap=2))

    if conflicts:
        body.append(dmc.Space(h=4))
        body.append(dmc.Stack(
            [dmc.Text(f"{c['channel']}: {c['counts']} → using {c['chosen']}",
                      size="xs", ff="monospace", c="dimmed")
             for c in conflicts],
            gap=2))

    if suggestions:
        body.append(dmc.Space(h=6))
        body.append(dmc.Text("Probably the same marker — tick to add an alias, "
                             "then load again", size="xs", fw=600))
        # Pattern-matching id: this group exists only once a report has been
        # rendered, and a callback keyed on a plain id would fail on every page
        # load before the first inspect.
        body.append(dmc.CheckboxGroup(
            id={"comp": "alias-suggestions"},
            value=[],
            children=dmc.Stack(
                [dmc.Checkbox(
                    label=f"{' = '.join(s['names'])} → {s['suggested']}",
                    # The value IS the CLI fragment, so nothing rebuilds it.
                    value=";".join(f"{n}={s['suggested']}"
                                   for n in s["names"] if n != s["suggested"]),
                    size="xs")
                 for s in suggestions],
                gap=2),
        ))
    return dmc.Stack(body, gap=4)


def panel(inspect: dict | None, excluded: list[str] | None = None):
    """The whole panel: three independently collapsible sections."""
    report = (inspect or {}).get("report")
    if not report:
        return dmc.Text("Load the workspaces to see what came in.",
                        size="sm", c="dimmed")

    excluded = list(excluded or [])
    workspaces = report.get("workspaces") or []
    # Problems first: someone opening this is looking for what went wrong, and
    # scrolling past twenty successes to find it is the failure mode of the log
    # this replaces.
    workspaces = sorted(workspaces,
                        key=lambda w: (w.get("status") == "loaded",
                                       w.get("label", "")))
    totals = report.get("totals") or {}
    markers = report.get("markers") or {}
    n_common = len(markers.get("common") or [])
    n_diverging = len(markers.get("not_common") or [])
    n_off = len(excluded)

    sections = []

    # --- workspaces, its own collapsible (the user asked for this) ----------
    ws_summary = f"{len(workspaces)} total"
    if n_off:
        ws_summary += f", {n_off} off"
    sections.append(dmc.AccordionItem(
        value="workspaces",
        children=[
            dmc.AccordionControl(dmc.Group(
                [dmc.Text("Workspaces", size="sm", fw=500),
                 dmc.Badge(ws_summary, size="xs", variant="light",
                           color="yellow" if n_off else "gray")],
                gap="xs")),
            dmc.AccordionPanel(dmc.Stack(
                [dmc.Group(
                    [dmc.Text("Switch one off to leave it out of the run.",
                              size="xs", c="dimmed"),
                     dmc.Button("Use all", id={"comp": "wsp-all"},
                                size="compact-xs", variant="subtle",
                                color="gray")],
                    justify="space-between")]
                + [_workspace_row(w, w["label"] in excluded)
                   for w in workspaces],
                gap=6)),
        ],
    ))

    # --- markers ------------------------------------------------------------
    sections.append(dmc.AccordionItem(
        value="markers",
        children=[
            dmc.AccordionControl(dmc.Group(
                [dmc.Text("Markers", size="sm", fw=500),
                 dmc.Badge(f"{n_common} shared" + (f", {n_diverging} not"
                                                   if n_diverging else ""),
                           size="xs", variant="light",
                           color="red" if not n_common else
                           ("yellow" if n_diverging else "gray"))],
                gap="xs")),
            dmc.AccordionPanel(marker_section(report)),
        ],
    ))

    # --- samples ------------------------------------------------------------
    unnamed = totals.get("samples_without_marker_names") or 0
    routes = totals.get("routes") or {}
    sample_lines = [
        dmc.Text(f"{totals.get('samples') or 0} of "
                 f"{totals.get('samples_listed') or 0} listed samples loaded.",
                 size="xs"),
    ]
    if routes.get("content"):
        sample_lines.append(dmc.Text(
            f"{routes['content']} paired by the $FIL keyword inside the file "
            f"rather than by name — they would be missing under strict "
            f"matching.", size="xs", c="dimmed"))
    if unnamed:
        sample_lines.append(dmc.Text(
            f"{unnamed} carry no marker names of their own; recovered from the "
            f"channel each was measured on.", size="xs", c="dimmed"))
    sections.append(dmc.AccordionItem(
        value="samples",
        children=[
            dmc.AccordionControl(dmc.Group(
                [dmc.Text("Samples", size="sm", fw=500),
                 dmc.Badge(f"{totals.get('samples') or 0}", size="xs",
                           variant="light", color="gray")],
                gap="xs")),
            dmc.AccordionPanel(dmc.Stack(sample_lines, gap=4)),
        ],
    ))

    return dmc.Accordion(
        multiple=True,
        # Nothing open by default. The headline badge already says whether
        # anything needs attention, and the callback in steps.py opens the
        # section that does.
        value=[],
        id="acc-health-sections",
        variant="separated",
        children=sections,
    )
