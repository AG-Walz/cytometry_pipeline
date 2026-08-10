"""Remember a workspace load, so reopening the app does not repeat it.

Reading 33 FlowJo workspaces builds a GatingSet per workspace and takes ten
minutes. Nothing about that load depends on the app being open, and nothing
about it changes between two launches pointed at the same folders — yet it was
repeated from scratch every time, which made the app feel broken in the only
place a user is already waiting.

This stores R's load report under a key derived from everything that could
change it, and hands it straight back when the key matches.

It computes NOTHING. The cached document is exactly the JSON
``src/visualization/umap_load_report.R`` wrote; this decides only whether it is
still the right answer.

WHAT THE KEY COVERS
    Every option that changes what the load produces — the folders, the sample
    group, the match mode, the marker aliases, the marker recovery switch — and
    the identity of the workspace FILES themselves (path, size, mtime). Edit a
    gate in FlowJo and the key changes, so the cache cannot serve a tree that no
    longer exists.

WHAT IT DELIBERATELY DOES NOT COVER
    The FCS files. There can be hundreds of thousands of them and stat-ing the
    tree would cost more than it saves. A cohort gains or loses an FCS file far
    less often than an app is reopened, and the Reload button exists for when
    it happens — which is why the button is not tucked away.
"""

from __future__ import annotations

import hashlib
import json
import os
import time
from pathlib import Path

from app import config

#: Where reports live. Separate from the form state so clearing one does not
#: silently clear the other.
CACHE_SUBDIR = "reports"

#: Reports older than this are ignored and re-read. A long-lived key can only
#: be wrong in ways the key does not cover (an FCS file appearing, a package
#: upgrade changing what CytoML accepts), and a week bounds how long a stale
#: answer can survive without anyone noticing.
MAX_AGE_SECONDS = 7 * 24 * 3600

#: Fields of the inspect request that change the report. A field NOT listed
#: here is one the cache will ignore — so anything added to the load must be
#: added here too, or the app will serve the previous answer for it.
#: The workspace SELECTION is deliberately not among them: a load reads every
#: workspace it finds (see rscript.inspect_argv), so flicking a switch in the
#: health panel cannot change the report and must not throw the cache away.
KEYED_FIELDS = (
    "wsp", "fcs_dir", "group", "match_mode", "marker_aliases",
    "recover_markers", "allow_missing_samples",
)


def _cache_dir() -> Path:
    return config.CACHE_DIR / CACHE_SUBDIR


def _workspace_fingerprint(wsp: str) -> list:
    """(path, size, mtime) for every .wsp under ``wsp``, sorted.

    Sorted so directory-order differences between filesystems cannot change the
    key. Unreadable entries contribute their error rather than being skipped:
    a workspace that has become unreadable is a different situation from one
    that is not there, and both must miss the cache.
    """
    path = Path(wsp)
    if path.is_file():
        files = [path]
    elif path.is_dir():
        files = sorted(path.rglob("*.wsp")) + sorted(path.rglob("*.WSP"))
    else:
        return [["missing", wsp]]

    out = []
    for f in sorted(set(files)):
        try:
            st = f.stat()
            out.append([str(f), st.st_size, int(st.st_mtime)])
        except OSError as exc:
            out.append([str(f), "unreadable", str(exc)])
    return out


def key_for(request: dict) -> str:
    """A stable hash of everything that decides what the load will report."""
    payload = {f: request.get(f) for f in KEYED_FIELDS}
    payload["_workspaces"] = _workspace_fingerprint(request.get("wsp") or "")
    blob = json.dumps(payload, sort_keys=True, default=str)
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()[:32]


def path_for(cache_key: str) -> Path:
    return _cache_dir() / f"{cache_key}.json"


def load(cache_key: str) -> dict | None:
    """The cached report, or None. Never raises: a bad cache is a cache miss."""
    p = path_for(cache_key)
    try:
        if not p.is_file():
            return None
        if time.time() - p.stat().st_mtime > MAX_AGE_SECONDS:
            return None
        return json.loads(p.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None


def reserve(cache_key: str) -> str:
    """The path R should write this load's report to.

    R writes straight into the cache rather than the app copying afterwards:
    one file, written once, and a run that dies half-way leaves a file that
    fails to parse and is therefore a miss — not a truncated report that looks
    valid.
    """
    d = _cache_dir()
    d.mkdir(parents=True, exist_ok=True)
    p = path_for(cache_key)
    try:
        p.unlink(missing_ok=True)
    except OSError:
        pass
    return str(p)


def age_text(cache_key: str) -> str:
    """How old the cached report is, for the UI to say so plainly."""
    try:
        seconds = time.time() - path_for(cache_key).stat().st_mtime
    except OSError:
        return ""
    if seconds < 90:
        return "just now"
    minutes = int(seconds // 60)
    if minutes < 60:
        return f"{minutes} min ago"
    hours = minutes // 60
    if hours < 24:
        return f"{hours} h ago"
    return f"{hours // 24} d ago"


def clear() -> int:
    """Drop every cached report. Returns how many were removed."""
    d = _cache_dir()
    if not d.is_dir():
        return 0
    n = 0
    for f in d.glob("*.json"):
        try:
            os.unlink(f)
            n += 1
        except OSError:
            pass
    return n
