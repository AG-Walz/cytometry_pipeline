"""Phase durations of complete runs, by settings: the reference for estimates.

One JSON file in the app's cache folder: {settings key: {"phases": {phase:
seconds}, "saved": epoch}}. The key is the run's argv without the output
folder and without --no_resume, because neither changes the work a run does:
the same settings into a new output folder take the same time.

A phase reused from a checkpoint took no time, so it is not written; the
value from an earlier run stays. Writes are best effort: an estimate is not
worth failing a run for.
"""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
import time
from pathlib import Path

from app import config

HISTORY_FILE = config.CACHE_DIR / "phase_history.json"

#: Arguments that do not change how long a run takes.
_FLAG_WITH_VALUE = {"--output_dir"}
_FLAG_ALONE = {"--no_resume"}

#: Entries kept. Old settings are dropped first.
MAX_ENTRIES = 200


def key_for(argv: list[str]) -> str:
    kept, skip = [], False
    for arg in argv:
        if skip:
            skip = False
            continue
        if arg in _FLAG_WITH_VALUE:
            skip = True
            continue
        if arg in _FLAG_ALONE or arg.split("=", 1)[0] in _FLAG_WITH_VALUE:
            continue
        kept.append(arg)
    return hashlib.sha1("\0".join(kept).encode("utf-8")).hexdigest()


def _read(path: Path) -> dict:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
        return data if isinstance(data, dict) else {}
    except (OSError, ValueError):
        return {}


def load(argv: list[str], path: Path | None = None) -> dict[str, float] | None:
    entry = _read(path or HISTORY_FILE).get(key_for(argv))
    if not isinstance(entry, dict):
        return None
    phases = entry.get("phases")
    if not isinstance(phases, dict) or not phases:
        return None
    try:
        return {str(k): float(v) for k, v in phases.items()}
    except (TypeError, ValueError):
        return None


def save(argv: list[str], phases: dict[str, float],
         path: Path | None = None) -> None:
    if not phases:
        return
    path = path or HISTORY_FILE
    data = _read(path)
    key = key_for(argv)
    old = data.get(key) if isinstance(data.get(key), dict) else {}
    merged = dict(old.get("phases") or {})
    merged.update({p: round(float(v), 1) for p, v in phases.items()})
    data[key] = {"phases": merged, "saved": time.time()}
    if len(data) > MAX_ENTRIES:
        newest = sorted(data.items(), key=lambda kv: kv[1].get("saved", 0),
                        reverse=True)[:MAX_ENTRIES]
        data = dict(newest)
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=".phase_history.")
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            json.dump(data, fh, indent=1)
        os.replace(tmp, path)
    except OSError:
        pass
