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

POLL_MS = 700


def _number(id_: str, label: str, value, **kwargs) -> dmc.NumberInput:
    return dmc.NumberInput(id=id_, label=label, value=value, **kwargs)


def _step_data(state: dict) -> list:
    return [
        dmc.Text(
            "Point at a FlowJo workspace and the directory holding its FCS "
            "files, then load it. The gating tree and marker panel only exist "
            "once the workspace has been read, so steps 2–4 stay locked until "
            "this succeeds.",
            size="sm",
            c="dimmed",
        ),
        dmc.Space(h=14),
        browser_field(
            "wsp",
            "FlowJo workspace (.wsp)",
            "The workspace whose gating you want to use",
            mode="file",
            extensions=config.WSP_EXTENSIONS,
            placeholder="/path/to/workspace.wsp",
            value=state["wsp"],
        ),
        dmc.Space(h=12),
        browser_field(
            "fcs",
            "FCS directory",
            "Directory containing the FCS files the workspace references",
            mode="dir",
            placeholder="/path/to/fcs",
            value=state["fcs_dir"],
        ),
        dmc.Space(h=12),
        dmc.Group(
            [
                _number(
                    "in-group",
                    "Workspace group",
                    state["group"],
                    min=1,
                    step=1,
                    allowDecimal=False,
                    description="Sample group index to load",
                    style={"flex": 1},
                ),
                _number(
                    "in-expected",
                    "Expected samples",
                    state["expected_samples"],
                    min=0,
                    step=1,
                    allowDecimal=False,
                    description="0 = no check; warns on mismatch",
                    style={"flex": 1},
                ),
            ],
            grow=True,
            align="flex-start",
        ),
        dmc.Space(h=16),
        dmc.Button(
            "Load workspace",
            id="btn-inspect",
            leftSection=dmc.Text("▸", size="sm"),
            fullWidth=True,
        ),
        dmc.Space(h=12),
        html.Div(id="inspect-feedback"),
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
    ]


def _step_params(state: dict) -> list:
    return [
        dmc.Text(
            "These map one-to-one onto the CLI options; anything left at its "
            "default is omitted from the command.",
            size="sm",
            c="dimmed",
        ),
        dmc.Space(h=14),
        dmc.Select(
            id="in-normalize",
            label="Per-marker normalisation",
            description="Computed from the pooled gated cells, before UMAP",
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
            description="0 = keep every cell; applied after gating",
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
        dmc.Space(h=12),
        dmc.Group(
            [
                dmc.Select(
                    id="in-metric",
                    label="metric",
                    data=[{"value": m, "label": m} for m in config.UMAP_METRICS],
                    value=state["metric"],
                    allowDeselect=False,
                    style={"flex": 1},
                ),
                _number("in-seed", "seed", state["seed"], min=0, step=1,
                        allowDecimal=False, style={"flex": 1}),
            ],
            grow=True,
            align="flex-start",
        ),
        dmc.Space(h=16),
        dmc.Divider(label="Plot appearance", labelPosition="center"),
        dmc.Space(h=10),
        dmc.Group(
            [
                _number("in-pointsize", "point_size", state["point_size"],
                        min=0.05, step=0.05, decimalScale=2, style={"flex": 1}),
                _number("in-pointalpha", "point_alpha", state["point_alpha"],
                        min=0.05, max=1, step=0.05, decimalScale=2,
                        style={"flex": 1}),
            ],
            grow=True,
            align="flex-start",
        ),
        dmc.Space(h=12),
        dmc.Group(
            [
                _number("in-width", "plot_width (in)", state["plot_width"],
                        min=1, step=1, decimalScale=1, style={"flex": 1}),
                _number("in-height", "plot_height (in)", state["plot_height"],
                        min=1, step=1, decimalScale=1, style={"flex": 1}),
            ],
            grow=True,
            align="flex-start",
        ),
        dmc.Space(h=16),
        dmc.CheckboxGroup(
            id="in-formats",
            label="Output formats",
            description=(
                "Exactly what is passed to --output_formats. PNG or SVG is "
                "needed for the in-app gallery; PDF-only runs show download "
                "links instead."
            ),
            value=state["output_formats"],
            children=dmc.Group(
                [dmc.Checkbox(label=f, value=f) for f in config.OUTPUT_FORMATS],
                mt=8,
            ),
        ),
        dmc.Space(h=8),
        html.Div(id="format-warning"),
    ]


def _step_run(state: dict) -> list:
    return [
        browser_field(
            "outparent",
            "Output parent directory",
            "The run writes into a new subfolder of this directory",
            mode="dir",
            placeholder="/path/to/results",
            value=state["output_parent"],
        ),
        dmc.Space(h=12),
        dmc.TextInput(
            id="in-outname",
            label="Output folder name",
            description="Pre-filled from the population and the current time",
            value=state["output_name"],
            placeholder="umap_viable_cells_20260728_1432",
        ),
        dmc.Space(h=8),
        html.Div(id="output-echo"),
        dmc.Space(h=16),
        dmc.Divider(label="Command", labelPosition="center"),
        dmc.Space(h=10),
        dmc.Text(
            "This is exactly what will be executed — copy it to reproduce the "
            "run from a shell.",
            size="xs",
            c="dimmed",
        ),
        dmc.Space(h=6),
        dmc.Code(id="cmd-preview", block=True, children="",
                 style={"whiteSpace": "pre-wrap", "wordBreak": "break-all"}),
        dmc.Space(h=14),
        html.Div(id="run-validation"),
        dmc.Space(h=10),
        dmc.Button("Run UMAP", id="btn-run", fullWidth=True, size="md"),
    ]


def _left_column(state: dict) -> dmc.Paper:
    return dmc.Paper(
        [
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
                                    # The per-marker accordion is a static
                                    # component so its open set can drive which
                                    # panels get built: Mantine keeps closed
                                    # panels mounted, so building them all
                                    # would put every marker x sample image on
                                    # the page at once.
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
                                dmc.Badge("umap_from_flowjo.R", variant="dot",
                                          color="gray", size="sm"),
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
                        dmc.GridCol(_left_column(state),
                                    span={"base": 12, "lg": 5}),
                        dmc.GridCol(_right_column(), span={"base": 12, "lg": 7}),
                    ],
                    gutter="md",
                ),
            ],
            p="md",
        ),
    )
