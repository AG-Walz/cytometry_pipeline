"""The page: a four-step stepper on the left, log and results on the right.

Every input component exists from the start (populated later by callbacks)
rather than being created on the fly. That keeps the callback graph static and
means a component can never be addressed before it exists.
"""

from __future__ import annotations

import dash_mantine_components as dmc
from dash import dcc, html

from app import config
from app.components.filebrowser import browser_field
from app.services import rscript

POLL_MS = 700


def _number(id_: str, label: str, value, **kwargs) -> dmc.NumberInput:
    return dmc.NumberInput(id=id_, label=label, value=value, **kwargs)


def _step_data(state: dict) -> list:
    """Step 1 — where the data is, and what came back.

    Two fields and a button. Everything that can be left alone is inside
    Advanced, and everything explaining WHY sits in a tooltip on the control it
    explains, rather than as a paragraph above it. The one exception is the
    Advanced badge: a setting that changes which samples load must never be
    able to hide behind a closed panel (see advanced_badge in callbacks/steps).
    """
    return [
        browser_field(
            "wsp",
            "FlowJo workspace folder",
            "Recursive; a single .wsp also works",
            mode="dir",
            extensions=config.WSP_EXTENSIONS,
            placeholder="/path/to/workspaces",
            value=state["wsp"],
        ),
        dmc.Space(h=6),
        html.Div(id="wsp-scan-preview"),
        dmc.Space(h=12),
        browser_field(
            "fcs",
            "FCS folder",
            "Recursive; any subfolder layout",
            mode="dir",
            placeholder="/path/to/fcs",
            value=state["fcs_dir"],
        ),
        dmc.Space(h=12),
        # Advanced stays ONE panel, deliberately: splitting it into four named
        # collapsibles replaces one closed header with four, which is more to
        # read and worry about, not less. It also stays in step 1 rather than
        # moving to Parameters, because the sample group decides which samples
        # load at all and must be set before the load. A closed Mantine
        # accordion keeps its panel mounted, so these ids stay addressable by
        # the form callbacks whether or not it is open.
        dmc.Accordion(
            id="acc-advanced",
            value=None,
            variant="contained",
            chevronPosition="left",
            children=[
                dmc.AccordionItem(
                    value="advanced",
                    children=[
                        dmc.AccordionControl(
                            dmc.Group(
                                [
                                    dmc.Text("Advanced", size="sm", fw=500),
                                    html.Div(id="advanced-badge"),
                                ],
                                gap="xs",
                                align="center",
                            )
                        ),
                        dmc.AccordionPanel(
                            [
                                dmc.TextInput(
                                    id="in-group",
                                    label="Sample group",
                                    value=str(state["group"]),
                                    description=(
                                        "Index or name. 1 is “All Samples” in "
                                        "most workspaces; prefer the name "
                                        "across several, since the index can "
                                        "differ. The Log lists each "
                                        "workspace's groups."
                                    ),
                                    placeholder="1",
                                ),
                                dmc.Space(h=12),
                                dmc.Switch(
                                    id="sw-allow-missing",
                                    label="Load anyway when FCS files are missing",
                                    description=(
                                        "Drops the samples whose .fcs file is "
                                        "missing instead of failing the run. "
                                        "The Log names every one."
                                    ),
                                    checked=state["allow_missing_samples"],
                                ),
                                dmc.Space(h=12),
                                dmc.RadioGroup(
                                    id="in-match-mode",
                                    value=state["match_mode"],
                                    label="Match samples to files",
                                    description=(
                                        "Renamed files no longer match by "
                                        "name; “by content” falls back to the "
                                        "$FIL keyword inside each file plus "
                                        "its event count, which a rename "
                                        "cannot change."
                                    ),
                                    children=dmc.Stack(
                                        [
                                            dmc.Radio(
                                                label="By content, then by name",
                                                value="content",
                                            ),
                                            dmc.Radio(
                                                label="By file name only",
                                                value="strict",
                                            ),
                                        ],
                                        gap=4,
                                    ),
                                    size="sm",
                                ),
                                dmc.Space(h=12),
                                # An input, so it lives with the inputs — not in
                                # the Data health panel that suggests values for
                                # it, which is rendered only after a load and
                                # would leave this unmounted on a fresh page.
                                dmc.TextInput(
                                    id="in-marker-aliases",
                                    label="Marker aliases",
                                    value=state["marker_aliases"],
                                    description=(
                                        "written=used, separated by “;”. "
                                        "Reconciles one marker spelled two "
                                        "ways. Data health suggests these."
                                    ),
                                    placeholder="INFy=IFN;TNFa=TNF",
                                ),
                                dmc.Space(h=12),
                                dmc.Switch(
                                    id="sw-recover-markers",
                                    label="Recover marker names from the channel",
                                    description=(
                                        "On: files with no $PnS names borrow "
                                        "them from files that have them, "
                                        "matched on detector."
                                    ),
                                    checked=state["recover_markers"],
                                ),
                            ]
                        ),
                    ],
                )
            ],
        ),
        dmc.Space(h=16),
        dmc.Group(
            [
                dmc.Button(
                    "Load workspaces",
                    id="btn-inspect",
                    leftSection=dmc.Text("▸", size="sm"),
                    style={"flex": 1},
                ),
                # The cache is keyed on the workspace files themselves, so it
                # cannot serve a tree that changed. It does NOT watch the FCS
                # folder — stat-ing hundreds of thousands of files would cost
                # more than it saves — so this is the way out when a file was
                # added there. Visible, not buried in Advanced, because it is
                # the answer to "why did nothing change?".
                dmc.Tooltip(
                    dmc.Button("Re-read", id="btn-reload", variant="subtle",
                               color="gray"),
                    label="Ignore the remembered result and read every "
                          "workspace again",
                    multiline=True, w=240, withArrow=True,
                ),
            ],
            gap="xs",
            grow=False,
        ),
        dmc.Space(h=12),
        html.Div(id="inspect-feedback"),
        dmc.Space(h=12),
        # Collapsed by default; the badge in the header carries the verdict, and
        # the callback in steps.py opens it only when something needs attention.
        dmc.Accordion(
            id="acc-health",
            value=None,
            variant="contained",
            chevronPosition="left",
            children=[
                dmc.AccordionItem(
                    value="health",
                    children=[
                        dmc.AccordionControl(
                            dmc.Group(
                                [
                                    dmc.Text("Data health", size="sm", fw=500),
                                    html.Div(id="health-badge"),
                                ],
                                gap="xs",
                            )
                        ),
                        dmc.AccordionPanel(html.Div(id="health-panel")),
                    ],
                )
            ],
        ),
    ]


def _step_gates(state: dict) -> list:
    return [
        html.Div(id="step2-summary"),
        dmc.Space(h=10),
        dmc.Group(
            [
                dmc.Text("Population", fw=600, size="sm"),
                dmc.Group(
                    [
                        dmc.Button("Expand all", id="btn-expand-all",
                                   size="compact-xs", variant="subtle"),
                        dmc.Button("Collapse", id="btn-collapse",
                                   size="compact-xs", variant="subtle"),
                    ],
                    gap=4,
                ),
            ],
            justify="space-between",
        ),
        dmc.Space(h=6),
        dmc.Switch(
            id="sw-per-wsp",
            label="Pick the population separately in each workspace",
            description=(
                "Use this when workspaces name the same biology differently "
                "— \"Live/CD8\" in one, \"Live cells/CD8\" in another. The "
                "list above them is only the gates ALL workspaces share, "
                "which collapses to a couple of nodes when they were gated "
                "independently. A workspace you leave unset is left out of "
                "the UMAP."
            ),
            checked=bool(state.get("population_map")),
        ),
        dmc.Space(h=8),
        # Always mounted, hidden when the mode is off: a callback whose Input
        # is missing from the layout never fires, and the gating tree depends
        # on this one.
        html.Div(
            id="per-wsp-select-wrap",
            style={"display": "none"},
            # value="" and never None: dash-mantine-components 2.8 fails to
            # render a Select whose value is None, and the error it raises is
            # a colour-parsing one, which points nowhere near the cause.
            children=dmc.Select(
                id="sel-wsp",
                label="Workspace",
                description=("The tree below shows this workspace's own gates; "
                             "clicking one sets its population."),
                data=[],
                value="",
                searchable=True,
                allowDeselect=False,
            ),
        ),
        dmc.Space(h=8),
        html.Div(id="per-wsp-picker"),
        html.Div(id="gate-tree-mode"),
        dmc.ScrollArea(
            dmc.Tree(
                id="gate-tree",
                data=[],
                expanded=[],
                selected=[],
                checked=[],
                selectOnClick=True,
                expandOnClick=True,
                checkboxes=False,
                allowRangeSelection=False,
            ),
            h=280,
            type="auto",
            offsetScrollbars=True,
            style={
                "border": "1px solid var(--mantine-color-default-border)",
                "borderRadius": "var(--mantine-radius-sm)",
                "padding": "8px",
            },
        ),
        dmc.Space(h=8),
        html.Div(id="population-echo"),
        dmc.Space(h=14),
        dmc.Divider(label="Subpopulation colouring", labelPosition="center"),
        dmc.Space(h=10),
        dmc.Switch(
            id="sw-explicit",
            label="Choose an explicit gate set instead of a mode",
            description=(
                "Tick gates in the tree above to build a custom colour "
                "partition (--subpop_select), overriding the mode below."
            ),
            checked=state["explicit_subpop"],
        ),
        dmc.Space(h=10),
        html.Div(
            id="subpop-mode-wrap",
            children=[
                dmc.SegmentedControl(
                    id="subpop-mode",
                    value=state["subpop_mode"],
                    data=[{"value": m, "label": m} for m in config.SUBPOP_MODES],
                    fullWidth=True,
                    size="xs",
                ),
                dmc.Space(h=6),
                html.Div(id="subpop-mode-help"),
            ],
        ),
        html.Div(id="subpop-chips"),
        dmc.Space(h=16),
        dmc.Divider(label="Markers", labelPosition="center"),
        dmc.Space(h=10),
        dmc.MultiSelect(
            id="in-markers",
            label="Markers to use",
            description=(
                "Leave empty to use every marker common to all samples. "
                "Selection order is preserved."
            ),
            data=[],
            value=state["markers"],
            searchable=True,
            clearable=True,
            hidePickedOptions=True,
            placeholder="all common markers",
            nothingFoundMessage="No such marker in this workspace",
        ),
        dmc.Space(h=16),
        dmc.Divider(label="Sample annotations", labelPosition="center"),
        dmc.Space(h=10),
        dmc.Text(
            "Group the samples by anything you know about them — patient, "
            "timepoint, response — and get one plot coloured by it plus one "
            "plot per level, all on the same embedding. Optional.",
            size="xs",
            c="dimmed",
        ),
        dmc.Space(h=10),
        browser_field(
            "meta",
            "Metadata CSV",
            "Keyed on sample; every other column becomes an annotation",
            mode="file",
            extensions=config.METADATA_EXTENSIONS,
            placeholder="/path/to/cohort.csv",
            value=state["sample_metadata"],
        ),
        dmc.Space(h=10),
        dmc.TextInput(
            id="in-pattern",
            label="…or derive from the sample name",
            description=(
                "A regex with named groups, used only when no CSV is given. "
                "Example: (?<patient>MTB[0-9]+)_(?<timepoint>TP[0-9]+)"
            ),
            value=state["sample_pattern"],
            placeholder="(?<patient>...)_(?<timepoint>...)",
        ),
        dmc.Space(h=8),
        html.Div(id="annotate-preview"),
        dmc.Space(h=10),
        dmc.MultiSelect(
            id="in-annotate",
            label="Annotations to plot",
            description="Leave empty to plot every column that resolves.",
            data=[],
            value=state["annotate"],
            searchable=True,
            clearable=True,
            hidePickedOptions=True,
            placeholder="all available",
            nothingFoundMessage="No such column",
        ),
    ]


def _step_params(state: dict) -> list:
    """Step 3 — the four knobs that change the result, and the rest behind one door.

    normalisation, downsample, n_neighbors and min_dist are real decisions:
    the first two change the embedding and the runtime, the other two are the
    UMAP parameters anyone has ever been told to tune. metric, seed and the four
    plot dimensions are defaults nobody changes — five of them only affect how
    the image looks and can be re-plotted, and seed only matters when you want
    to reproduce, at which point the command preview in step 4 is what you copy.
    """
    return [
        dmc.Select(
            id="in-normalize",
            label="Per-marker normalisation",
            description="From the pooled gated cells, before UMAP",
            data=[{"value": n, "label": n} for n in config.NORMALIZE_METHODS],
            value=state["normalize"],
            allowDeselect=False,
        ),
        dmc.Space(h=12),
        _number(
            "in-downsample",
            "Downsample (cells per sample)",
            state["downsample"],
            min=0,
            step=1000,
            allowDecimal=False,
            description="0 = every cell",
        ),
        dmc.Space(h=12),
        dmc.Group(
            [
                _number("in-neighbors", "n_neighbors", state["n_neighbors"],
                        min=2, step=1, allowDecimal=False, style={"flex": 1}),
                _number("in-mindist", "min_dist", state["min_dist"],
                        min=0, max=1, step=0.05, decimalScale=3,
                        style={"flex": 1}),
            ],
            grow=True,
            align="flex-start",
        ),
        dmc.Space(h=16),
        dmc.CheckboxGroup(
            id="in-formats",
            label="Output formats",
            description="PNG or SVG for the in-app gallery.",
            value=state["output_formats"],
            children=dmc.Group(
                [dmc.Checkbox(label=f, value=f) for f in config.OUTPUT_FORMATS],
                mt=8,
            ),
        ),
        dmc.Space(h=8),
        # Stays visible: it reports that the run would produce no files.
        html.Div(id="format-warning"),
        dmc.Space(h=14),
        # A closed Mantine panel stays mounted, so every id below is still
        # addressable by the form callbacks. The badge is not decoration — it is
        # the rule from advanced_badge: a non-default setting may not hide.
        dmc.Accordion(
            id="acc-params",
            value=None,
            variant="contained",
            chevronPosition="left",
            children=[
                dmc.AccordionItem(
                    value="appearance",
                    children=[
                        dmc.AccordionControl(
                            dmc.Group(
                                [
                                    dmc.Text("Appearance & advanced", size="sm",
                                             fw=500),
                                    html.Div(id="params-badge"),
                                ],
                                gap="xs",
                            )
                        ),
                        dmc.AccordionPanel(
                            [
                                dmc.Group(
                                    [
                                        dmc.Select(
                                            id="in-metric",
                                            label="metric",
                                            data=[{"value": m, "label": m}
                                                  for m in config.UMAP_METRICS],
                                            value=state["metric"],
                                            allowDeselect=False,
                                            style={"flex": 1},
                                        ),
                                        _number("in-seed", "seed", state["seed"],
                                                min=0, step=1,
                                                allowDecimal=False,
                                                style={"flex": 1}),
                                    ],
                                    grow=True, align="flex-start",
                                ),
                                dmc.Space(h=12),
                                dmc.Group(
                                    [
                                        _number("in-pointsize", "point_size",
                                                state["point_size"], min=0.05,
                                                step=0.05, decimalScale=2,
                                                style={"flex": 1}),
                                        _number("in-pointalpha", "point_alpha",
                                                state["point_alpha"], min=0.05,
                                                max=1, step=0.05,
                                                decimalScale=2,
                                                style={"flex": 1}),
                                    ],
                                    grow=True, align="flex-start",
                                ),
                                dmc.Space(h=12),
                                dmc.Group(
                                    [
                                        _number("in-width", "plot_width (in)",
                                                state["plot_width"], min=1,
                                                step=1, decimalScale=1,
                                                style={"flex": 1}),
                                        _number("in-height", "plot_height (in)",
                                                state["plot_height"], min=1,
                                                step=1, decimalScale=1,
                                                style={"flex": 1}),
                                    ],
                                    grow=True, align="flex-start",
                                ),
                            ]
                        ),
                    ],
                )
            ],
        ),
    ]


def _step_run(state: dict) -> list:
    """Step 4 — where it goes, and the button.

    The command preview moves behind a disclosure. It is the reproducibility
    escape hatch, not something anyone reads before every run — and it is long
    enough to push the Run button below the fold. What stays visible is the
    output echo (which warns that a directory already exists and files will be
    overwritten) and the validation list.
    """
    return [
        browser_field(
            "outparent",
            "Output parent directory",
            "A new subfolder is created here",
            mode="dir",
            placeholder="/path/to/results",
            value=state["output_parent"],
        ),
        dmc.Space(h=12),
        dmc.TextInput(
            id="in-outname",
            label="Output folder name",
            value=state["output_name"],
            placeholder="umap_viable_cells_20260728_1432",
        ),
        dmc.Space(h=8),
        # Stays visible: this is where "the directory exists and files will be
        # overwritten" is said, and it is the only destructive warning here.
        html.Div(id="output-echo"),
        dmc.Space(h=14),
        dmc.Accordion(
            id="acc-command",
            value=None,
            variant="contained",
            chevronPosition="left",
            children=[
                dmc.AccordionItem(
                    value="command",
                    children=[
                        dmc.AccordionControl(
                            dmc.Text("Command", size="sm", fw=500)),
                        dmc.AccordionPanel(
                            dmc.Code(id="cmd-preview", block=True, children="",
                                     style={"whiteSpace": "pre-wrap",
                                            "wordBreak": "break-all"})),
                    ],
                )
            ],
        ),
        dmc.Space(h=14),
        html.Div(id="run-validation"),
        dmc.Space(h=10),
        dmc.Button("Run UMAP", id="btn-run", fullWidth=True, size="md"),
    ]


def _fcs_names_panel(state: dict) -> list:
    """The second tool: make each FCS file's $FIL match its name on disk.

    Controls only. The report streams into the shared Log pane on the right,
    which is why this tab bar wraps the left column and not the whole grid.

    The mode picker is UNCHANGED — same five modes, same read-only-first order,
    same default. That ordering is the safety design of this tab (app/README.md)
    and nothing about it is hidden, shortened or reordered here; only the prose
    around it is.
    """
    return [
        dmc.Badge("sync_fcs_names.R", variant="dot", color="gray", size="sm"),
        dmc.Space(h=10),
        dmc.Text(
            "FlowJo shows each file's $FIL keyword, not its name on disk. "
            "This writes the name on disk into $FIL.",
            size="sm", c="dimmed",
        ),
        dmc.Space(h=8),
        dmc.Accordion(
            id="acc-fcs-help",
            value=None,
            variant="contained",
            chevronPosition="left",
            children=[
                dmc.AccordionItem(
                    value="why",
                    children=[
                        dmc.AccordionControl(
                            dmc.Text("Why the names look wrong", size="sm")),
                        dmc.AccordionPanel(
                            dmc.Text(
                                "A run acquired under one file name records "
                                "that name inside every file, in $FIL. "
                                "Renaming the files on disk does not change "
                                "it, so the workspace shows a column of "
                                "identical samples. Check FlowJo → Preferences "
                                "→ Workspace → “How should data files be "
                                "named” first: set to “Use Short File Name”, "
                                "you do not need this at all.",
                                size="xs", c="dimmed")),
                    ],
                )
            ],
        ),
        dmc.Space(h=14),
        browser_field(
            "fcssync",
            "FCS folder",
            "Searched recursively",
            mode="dir",
            placeholder="/path/to/renamed",
            value=state.get("fcs_sync_dir", ""),
        ),
        dmc.Space(h=14),
        dmc.RadioGroup(
            id="fcs-mode",
            value="probe",
            label="What to do",
            # The escalation, attached to the control that orders it rather
            # than left in a paragraph above.
            description="Census → dry run → apply.",
            children=dmc.Stack(
                [
                    dmc.Radio(value=k, label=v["label"])
                    for k, v in rscript.SYNC_MODES.items()
                ],
                gap=6,
                mt=6,
            ),
        ),
        dmc.Space(h=12),
        # Always mounted, never conditionally rendered: hiding it would need a
        # callback, and layout.py's doctrine is that every input exists from
        # the start. It is simply ignored by every mode but "copy".
        browser_field(
            "fcsout",
            "Output folder",
            "Copies mode only, outside the FCS folder",
            mode="dir",
            placeholder="/path/to/patched_copies",
        ),
        dmc.Space(h=12),
        dmc.Accordion(
            id="acc-fcs-more",
            value=None,
            variant="contained",
            chevronPosition="left",
            children=[
                dmc.AccordionItem(
                    value="more",
                    children=[
                        dmc.AccordionControl(
                            dmc.Group(
                                [dmc.Text("More options", size="sm"),
                                 html.Div(id="fcs-also-badge")],
                                gap="xs")),
                        dmc.AccordionPanel(
                            dmc.TextInput(
                                id="fcs-also",
                                label="Also sync these keywords",
                                description=(
                                    "Comma-separated, e.g. $SRC,TUBE NAME. "
                                    "$FIL is always written; the replaced name "
                                    "is kept in ORIGINALFIL."
                                ),
                                placeholder="",
                            )),
                    ],
                )
            ],
        ),
        dmc.Space(h=16),
        # Visible, always: this is the pre-write disclosure, and the documented
        # way to reach the flags that have no UI on purpose.
        dmc.Divider(label="Command", labelPosition="center"),
        dmc.Space(h=8),
        dmc.Code(id="fcs-cmd", block=True, children="",
                 style={"whiteSpace": "pre-wrap", "wordBreak": "break-all"}),
        dmc.Space(h=10),
        html.Div(id="fcs-note"),
        dmc.Space(h=8),
        dmc.Button("Run census", id="btn-fcs", fullWidth=True, size="md"),
        dmc.Space(h=8),
        html.Div(id="fcs-feedback"),
        dmc.Space(h=18),
        dmc.Divider(label="Stored metadata", labelPosition="center"),
        dmc.Space(h=8),
        dmc.Group(
            [
                dmc.Tooltip(
                    dmc.Button("Show metadata", id="btn-fcs-kw",
                               variant="default", size="compact-sm"),
                    label="Decodes every keyword of every file in the folder "
                          "so the cohort can be compared. Read-only.",
                    multiline=True, w=250, withArrow=True,
                ),
                dmc.Switch(id="sw-fcs-allkw", label="Every keyword",
                           checked=False, size="xs"),
            ],
            justify="space-between", align="center",
        ),
        dmc.Space(h=10),
        dcc.Store(id="store-fcs-kw", data=None),
        html.Div(id="fcs-kw-table"),
    ]
def _left_column(state: dict) -> dmc.Paper:
    return dmc.Paper(
        [
            dmc.Badge("umap_from_flowjo.R", variant="dot", color="gray",
                      size="sm", mb="sm"),
            dmc.Stepper(
                id="stepper",
                active=0,
                allowNextStepsSelect=False,
                size="sm",
                iconSize=32,
                # Labels only: with descriptions the four steps wrap onto two
                # rows in this column width. Each pane opens with its own
                # explanatory line instead.
                children=[
                    dmc.StepperStep(
                        label="Data",
                        children=dmc.Box(_step_data(state), mt="md"),
                    ),
                    dmc.StepperStep(
                        label="Population",
                        children=dmc.Box(_step_gates(state), mt="md"),
                    ),
                    dmc.StepperStep(
                        label="Parameters",
                        children=dmc.Box(_step_params(state), mt="md"),
                    ),
                    dmc.StepperStep(
                        label="Run",
                        children=dmc.Box(_step_run(state), mt="md"),
                    ),
                ],
            ),
            dmc.Space(h=18),
            dmc.Group(
                [
                    dmc.Button("Back", id="btn-back", variant="default"),
                    dmc.Button("Next", id="btn-next"),
                ],
                justify="space-between",
            ),
        ],
        p="lg",
        withBorder=True,
        radius="md",
    )


def _right_column() -> dmc.Paper:
    return dmc.Paper(
        [
            dmc.Tabs(
                [
                    dmc.TabsList(
                        [
                            dmc.TabsTab("Log", value="log"),
                            dmc.TabsTab("Results", value="results"),
                        ]
                    ),
                    dmc.TabsPanel(
                        [
                            dmc.Space(h=12),
                            dmc.Group(
                                [
                                    dmc.Group(
                                        [
                                            dmc.Badge("idle", id="status-badge",
                                                      color="gray", variant="light"),
                                            dmc.Text("—", id="elapsed-text",
                                                     size="sm", c="dimmed"),
                                        ],
                                        gap="xs",
                                    ),
                                    dmc.Group(
                                        [
                                            dmc.Switch(
                                                id="sw-follow",
                                                label="Follow",
                                                checked=True,
                                                size="xs",
                                            ),
                                            dmc.Button(
                                                "Cancel",
                                                id="btn-cancel",
                                                color="red",
                                                variant="light",
                                                size="compact-sm",
                                                disabled=True,
                                            ),
                                        ],
                                        gap="sm",
                                    ),
                                ],
                                justify="space-between",
                            ),
                            # Progress. Empty until R reports a phase, and
                            # hidden again when nothing is running — a bar
                            # frozen at 40% after a finished run reads as a
                            # stalled job.
                            html.Div(id="progress-wrap"),
                            dmc.Space(h=10),
                            dmc.ScrollArea(
                                html.Pre(
                                    "Nothing has run yet.",
                                    id="log-text",
                                    style={
                                        "margin": 0,
                                        "fontSize": "12px",
                                        "lineHeight": "1.45",
                                        "whiteSpace": "pre-wrap",
                                        "wordBreak": "break-word",
                                    },
                                ),
                                id="log-scroll",
                                h="calc(100vh - 260px)",
                                type="auto",
                                offsetScrollbars=True,
                                style={
                                    "border": "1px solid var(--mantine-color-default-border)",
                                    "borderRadius": "var(--mantine-radius-sm)",
                                    "padding": "10px",
                                    "background": "var(--mantine-color-body)",
                                },
                            ),
                        ],
                        value="log",
                    ),
                    dmc.TabsPanel(
                        [
                            dmc.Space(h=12),
                            dmc.Group(
                                [
                                    dmc.Text("No results yet.", id="results-path",
                                             size="sm", c="dimmed",
                                             style={"wordBreak": "break-all",
                                                    "flex": 1}),
                                    dmc.Button("Refresh", id="btn-refresh-results",
                                               variant="default",
                                               size="compact-sm"),
                                ],
                                justify="space-between",
                                align="center",
                                wrap="nowrap",
                            ),
                            dmc.Space(h=8),
                            dmc.Group(
                                [
                                    dmc.Select(
                                        id="filter-marker",
                                        data=[],
                                        placeholder="all markers",
                                        label="Marker",
                                        clearable=True,
                                        searchable=True,
                                        size="xs",
                                        style={"flex": 1},
                                    ),
                                    dmc.Select(
                                        id="filter-sample",
                                        data=[],
                                        placeholder="all samples",
                                        label="Sample",
                                        clearable=True,
                                        searchable=True,
                                        size="xs",
                                        style={"flex": 1},
                                    ),
                                ],
                                grow=True,
                                align="flex-end",
                            ),
                            dmc.Space(h=10),
                            dmc.ScrollArea(
                                [
                                    html.Div(id="gallery"),
                                    # Both accordions are static components so
                                    # their open sets can drive which panels
                                    # get built: Mantine keeps closed panels
                                    # mounted, so building them all would put
                                    # every image on the page at once — marker
                                    # x sample, and every annotation level.
                                    dmc.Accordion(
                                        id="annot-accordion",
                                        multiple=True,
                                        value=[],
                                        variant="separated",
                                        children=[],
                                        mb="md",
                                    ),
                                    html.Div(id="per-marker-header"),
                                    dmc.Accordion(
                                        id="marker-accordion",
                                        multiple=True,
                                        value=[],
                                        variant="separated",
                                        children=[],
                                    ),
                                ],
                                h="calc(100vh - 300px)",
                                type="auto",
                                offsetScrollbars=True,
                            ),
                        ],
                        value="results",
                    ),
                ],
                id="tabs",
                value="log",
            ),
        ],
        p="lg",
        withBorder=True,
        radius="md",
    )


def build_layout(state: dict):
    return dmc.MantineProvider(
        id="mantine-provider",
        forceColorScheme="light",
        children=dmc.Box(
            [
                dcc.Store(id="store-inspect", data={}),
                # {workspace: gate}. A Store rather than a component prop
                # because the mapping is per workspace and there is no single
                # widget that owns it; it is read back through FORM_FIELDS like
                # any other field, and persists with the rest of the form.
                dcc.Store(id="store-popmap", data=state.get("population_map") or {}),
                # Workspace labels the user has switched OFF. An exclusion list
                # rather than an inclusion one so the default — everything that
                # loaded — needs no state at all, and so a workspace appearing
                # in the folder later is included rather than silently dropped.
                dcc.Store(id="store-wsp-exclude",
                          data=state.get("excluded_workspaces") or []),
                dcc.Store(id="store-results-dir", data=""),
                dcc.Store(id="store-run-generation", data=-1),
                dcc.Interval(id="poll", interval=POLL_MS, disabled=False),
                # Sink for callbacks whose only real effect is starting or
                # stopping the subprocess; Dash requires every callback to have
                # an output.
                html.Div(id="run-start-sink", style={"display": "none"}),
                dmc.Modal(id="lightbox", size="90%", opened=False, zIndex=3000,
                          children=html.Div(id="lightbox-body")),
                dmc.Group(
                    [
                        dmc.Group(
                            [
                                dmc.Title("FlowJo → UMAP", order=3),
                            ],
                            gap="sm",
                            align="center",
                        ),
                        dmc.Group(
                            [
                                dmc.Button("Reset form", id="btn-reset",
                                           variant="subtle", color="gray",
                                           size="compact-sm"),
                                dmc.Switch(id="sw-theme", label="Dark",
                                           checked=False, size="sm"),
                            ],
                            gap="md",
                        ),
                    ],
                    justify="space-between",
                    align="center",
                    mb="md",
                ),
                dmc.Grid(
                    [
                        # The tab bar wraps the LEFT column only, so the Log
                        # and Results panes stay visible for both tools. Tabs
                        # around the whole Grid would hide the log behind the
                        # UMAP tab — and the log is where the FCS tool's entire
                        # output goes. id is "tool-tabs": "tabs" belongs to the
                        # right column and is both a State and an Output of the
                        # poll callback.
                        dmc.GridCol(
                            dmc.Tabs(
                                [
                                    dmc.TabsList([
                                        dmc.TabsTab("UMAP", value="umap"),
                                        dmc.TabsTab("FCS names", value="fcsnames"),
                                    ]),
                                    dmc.TabsPanel(_left_column(state),
                                                  value="umap", keepMounted=True),
                                    dmc.TabsPanel(
                                        dmc.Paper(_fcs_names_panel(state), p="lg",
                                                  withBorder=True, radius="md",
                                                  mt="md"),
                                        value="fcsnames", keepMounted=True),
                                ],
                                id="tool-tabs", value="umap",
                            ),
                            span={"base": 12, "lg": 5}),
                        dmc.GridCol(_right_column(), span={"base": 12, "lg": 7}),
                    ],
                    gutter="md",
                ),
            ],
            p="md",
        ),
    )
