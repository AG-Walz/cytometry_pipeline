"""Filesystem listing helpers for the dmc.Tree file/directory browser.

The tree is rooted at ``/`` and children are listed lazily: an unexpanded
directory carries a single placeholder child so Mantine draws its chevron, and
the placeholder is swapped for the real listing the first time the node is
expanded. That keeps the payload small on directories with many entries.
"""

from __future__ import annotations

import os
from pathlib import Path

#: Prefixes marking a not-yet-loaded / non-selectable node. Every real node
#: value is an absolute path, so anything not starting with "/" is a sentinel —
#: a test that survives JSON round-tripping and is trivial to mirror in JS.
PLACEHOLDER_PREFIX = "::placeholder::"
ERROR_PREFIX = "::error::"
EMPTY_PREFIX = "::empty::"
#: A file shown for orientation inside a *directory* picker. You are choosing a
#: folder, but seeing the .wsp files in it is how you know it is the right one,
#: so they are listed and not selectable.
INFO_PREFIX = "::info::"


def is_sentinel(value: str | None) -> bool:
    """True for placeholder / error / empty marker nodes (never selectable)."""
    return not value or not value.startswith("/")


def _placeholder(path: str) -> dict:
    return {"value": PLACEHOLDER_PREFIX + path, "label": "loading…"}


def _dir_node(path: Path) -> dict:
    return {
        "value": str(path),
        "label": path.name or str(path),
        "children": [_placeholder(str(path))],
    }


def _file_node(path: Path) -> dict:
    return {"value": str(path), "label": path.name}


def list_children(
    path: str,
    dirs_only: bool = False,
    extensions: tuple[str, ...] | None = None,
    show_hidden: bool = False,
) -> list[dict]:
    """Return dmc.Tree child nodes for ``path``.

    Directories always come first, then files, each sorted case-insensitively.
    Unreadable directories yield a single non-selectable error node rather than
    raising, so one bad directory cannot break the whole browser.
    """
    p = Path(path)
    try:
        entries = list(os.scandir(p))
    except (PermissionError, OSError) as exc:
        return [{"value": ERROR_PREFIX + str(p), "label": f"⚠ {exc.strerror or exc}"}]

    dirs: list[Path] = []
    files: list[Path] = []
    for entry in entries:
        if not show_hidden and entry.name.startswith("."):
            continue
        try:
            # follow_symlinks=True so a symlinked data directory still expands.
            if entry.is_dir():
                dirs.append(Path(entry.path))
            elif entry.is_file():
                files.append(Path(entry.path))
        except OSError:
            # Broken symlink or a race with a deleted entry — skip it.
            continue

    dirs.sort(key=lambda x: x.name.lower())
    files.sort(key=lambda x: x.name.lower())

    if extensions:
        wanted = tuple(e.lower() for e in extensions)
        files = [f for f in files if f.suffix.lower() in wanted]

    nodes = [_dir_node(d) for d in dirs]
    if not dirs_only:
        nodes.extend(_file_node(f) for f in files)
    elif extensions:
        # Directory mode with a filter: show the matching files so you can tell
        # you are standing in the right folder, but keep them unselectable —
        # is_sentinel() and validate_selection() both key on the "::" prefix.
        nodes.extend(
            {"value": INFO_PREFIX + str(f), "label": f.name} for f in files
        )

    if not nodes:
        label = "(no subdirectories)" if dirs_only else "(nothing matching here)"
        return [{"value": EMPTY_PREFIX + str(p), "label": label}]
    return nodes


def find_files(
    root: str, extensions: tuple[str, ...], limit: int = 500
) -> tuple[list[str], bool]:
    """Every file under ``root`` matching ``extensions``, searched recursively.

    Mirrors what the R side does with ``--wsp`` and ``--fcs_dir`` so the form
    can show you what a folder contains *before* you spend a workspace load
    finding out. Returns ``(paths, truncated)``; ``root`` may also be a single
    matching file, which comes back as a one-element list.

    Unreadable subdirectories are skipped rather than raising — the browser is
    rooted at ``/`` and will meet plenty of them.
    """
    p = Path(root)
    wanted = tuple(e.lower() for e in extensions)
    if p.is_file():
        return ([str(p)] if p.suffix.lower() in wanted else [], False)
    if not p.is_dir():
        return [], False

    found: list[str] = []
    for dirpath, dirnames, filenames in os.walk(p, onerror=lambda _e: None):
        dirnames[:] = [d for d in dirnames if not d.startswith(".")]
        for name in filenames:
            if name.lower().endswith(wanted):
                found.append(os.path.join(dirpath, name))
                if len(found) > limit:
                    return sorted(found[:limit]), True
    return sorted(found), False


def root_data() -> list[dict]:
    """Initial tree data: a single unexpanded ``/`` node."""
    return [{"value": "/", "label": "/", "children": [_placeholder("/")]}]


def find_node(nodes: list[dict], value: str) -> dict | None:
    """Depth-first search for the node carrying ``value``."""
    for node in nodes:
        if node.get("value") == value:
            return node
        children = node.get("children")
        if children:
            hit = find_node(children, value)
            if hit is not None:
                return hit
    return None


def node_needs_loading(node: dict | None) -> bool:
    """True when ``node`` still holds only its placeholder child."""
    if not node:
        return False
    children = node.get("children")
    if not children:
        return False
    return len(children) == 1 and str(children[0].get("value", "")).startswith(
        PLACEHOLDER_PREFIX
    )


def load_into(
    nodes: list[dict],
    value: str,
    dirs_only: bool,
    extensions: tuple[str, ...] | None,
    show_hidden: bool = False,
) -> bool:
    """Populate the children of ``value`` in place. Returns True if it changed."""
    node = find_node(nodes, value)
    if not node_needs_loading(node):
        return False
    node["children"] = list_children(value, dirs_only, extensions, show_hidden)
    return True


def ancestors(path: str) -> list[str]:
    """Every ancestor of ``path`` from ``/`` downwards, excluding ``path``."""
    p = Path(path)
    chain = [str(a) for a in reversed(p.parents)]
    return chain


def expand_to(
    path: str,
    dirs_only: bool,
    extensions: tuple[str, ...] | None,
    show_hidden: bool = False,
) -> tuple[list[dict], list[str]]:
    """Build tree data with every level down to ``path`` already loaded.

    Returns ``(data, expanded)``. Used when opening the browser on a field that
    already holds a path, and by the "jump to path" box, so the user lands on
    their directory instead of at ``/``. A path that does not exist degrades to
    the deepest existing ancestor.
    """
    data = root_data()
    target = Path(path)
    # Walk down to the deepest ancestor that actually exists.
    while str(target) != target.root and not target.exists():
        target = target.parent

    chain = ancestors(str(target)) + [str(target)]
    expanded: list[str] = []
    for value in chain:
        if not Path(value).is_dir():
            break
        load_into(data, value, dirs_only, extensions, show_hidden)
        expanded.append(value)
    return data, expanded


def validate_selection(value: str, dirs_only: bool, extensions: tuple[str, ...] | None):
    """Check a tree selection against the browser's mode.

    Returns ``(ok, message)`` — ``message`` explains the rejection, or describes
    the accepted selection.
    """
    if value and value.startswith(INFO_PREFIX):
        return False, "That is a file — select the folder that contains it."
    if not value or is_sentinel(value):
        return False, "Select a " + ("directory" if dirs_only else "file") + "."
    p = Path(value)
    if dirs_only:
        if not p.is_dir():
            return False, "That is a file — select a directory."
        return True, str(p)
    if p.is_dir():
        return False, "That is a directory — select a file."
    if extensions and p.suffix.lower() not in tuple(e.lower() for e in extensions):
        wanted = ", ".join(extensions)
        return False, f"Select a {wanted} file."
    return True, str(p)
