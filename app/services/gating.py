"""Turn the flat population paths from the workspace into dmc.Tree data.

``gs_get_pop_paths(path = "full")`` returns absolute gate paths such as
``/Lymphocytes/Single Cells/Viable cells``. The tree nests them on ``/`` and
keeps the full path as each node's value, so a selection is directly usable as
``--population`` and a checked node is directly usable for ``--subpop_select``.
"""

from __future__ import annotations


def build_tree(paths: list[str], labels: dict | None = None) -> list[dict]:
    """Nested dmc.Tree data covering every path, parents included.

    A path with no "/" — in practice the bare ``root`` node the gating tree
    starts with — is kept verbatim as its own value. Splitting it like the rest
    would invent "/root", which is not a gate and would not resolve on the R
    side.

    ``labels`` maps a full gate path to the text to show for it, and comes from
    R (``gate_display_label``). FlowJo permits a gate named with nothing but
    spaces, and real workspaces contain them — two siblings called " " and "  "
    render as two invisible, indistinguishable rows. R supplies a readable
    stand-in; the VALUE stays the untouched path, because that is what
    ``--population`` has to resolve against.

    Splitting on "/" keeps empty segments out, so an intermediate node invented
    for a parent path still falls back to its own segment text.
    """
    labels = labels or {}
    root: dict = {}
    for path in paths:
        if "/" not in path:
            root.setdefault(path, {"__value__": path, "__children__": {}})
            continue
        parts = [p for p in path.split("/") if p]
        node = root
        for i, part in enumerate(parts):
            full = "/" + "/".join(parts[: i + 1])
            node = node.setdefault(part, {"__value__": full, "__children__": {}})[
                "__children__"
            ]

    def to_nodes(mapping: dict) -> list[dict]:
        out = []
        for segment, payload in mapping.items():
            value = payload["__value__"]
            node: dict = {"value": value,
                          "label": labels.get(value) or segment}
            children = to_nodes(payload["__children__"])
            if children:
                node["children"] = children
            out.append(node)
        return out

    return to_nodes(root)


def default_expanded(paths: list[str], depth: int = 2) -> list[str]:
    """Ancestor paths up to ``depth`` levels, so the tree opens partly unfolded."""
    expanded: set[str] = set()
    for path in paths:
        if "/" not in path:      # the bare "root" node has nothing to unfold
            continue
        parts = [p for p in path.split("/") if p]
        for i in range(min(depth, len(parts))):
            expanded.add("/" + "/".join(parts[: i + 1]))
    return sorted(expanded)


def all_values(paths: list[str]) -> list[str]:
    """Every node value in the tree, including intermediate parents."""
    values: set[str] = set()
    for path in paths:
        if "/" not in path:      # "root" is its own value, not "/root"
            values.add(path)
            continue
        parts = [p for p in path.split("/") if p]
        for i in range(len(parts)):
            values.add("/" + "/".join(parts[: i + 1]))
    return sorted(values)


def descendants(paths: list[str], population: str) -> list[str]:
    """Strict descendants of ``population`` among ``paths``."""
    if not population:
        return []
    prefix = population.rstrip("/") + "/"
    return [p for p in all_values(paths) if p.startswith(prefix)]


def relative(path: str, population: str) -> str:
    """``path`` expressed relative to ``population`` (as --subpop_select wants)."""
    prefix = population.rstrip("/") + "/"
    return path[len(prefix):] if path.startswith(prefix) else path


def flat_subgate_tree(paths: list[str], population: str) -> list[dict]:
    """A FLAT, checkable list of every gate under ``population``.

    Mantine's Tree reports only *leaf* nodes in ``checked``: ticking a parent
    yields its children instead of the parent itself, so in a nested tree an
    intermediate gate such as ``CD4+`` could never be chosen as a colour class.
    Flattening makes every gate a leaf, so ``checked`` is exactly what the user
    ticked. Labels keep the relative path, so the hierarchy stays readable.
    """
    return [
        {"value": path, "label": relative(path, population)}
        for path in descendants(paths, population)
    ]
