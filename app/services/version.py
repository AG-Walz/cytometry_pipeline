"""Which code this app is running, and whether newer code exists.

The app runs whatever was checked out when its process started, and keeps
running it: a `git pull` on the server changes the files, not the process. On
2026-11-03 an instance started at 12:02 kept serving the old form for hours
after five fixes had been merged and pulled, and nothing on the page said so.

Three commits are compared:

    running  the checkout when this process started (read once, at import)
    on disk  the checkout now          -> differs: restart the app to load it
    GitHub   origin's main branch now  -> differs: pull, then restart

GitHub is asked at most every CHECK_SECONDS and with a short timeout; when it
cannot be reached the badge says so instead of guessing.
"""

from __future__ import annotations

import subprocess
import threading
import time

from app.config import REPO_ROOT

#: How often GitHub is asked for its main branch, at most.
CHECK_SECONDS = 600
#: How long one git call may take before it is given up on.
TIMEOUT_SECONDS = 8


def _git(*args: str, timeout: float = TIMEOUT_SECONDS) -> str | None:
    try:
        out = subprocess.run(["git", "-C", str(REPO_ROOT), *args],
                             capture_output=True, text=True, timeout=timeout,
                             check=True)
    except (OSError, subprocess.SubprocessError):
        return None
    return out.stdout.strip()


def _describe(ref: str = "HEAD") -> dict:
    """commit, short, date (YYYY-MM-DD HH:MM) and subject of ``ref``."""
    line = _git("log", "-1", "--format=%H%x1f%h%x1f%ci%x1f%s", ref)
    if not line:
        return {}
    commit, short, date, subject = (line.split("\x1f") + [""] * 4)[:4]
    return {"commit": commit, "short": short, "date": date[:16],
            "subject": subject}


#: The code this process loaded. Read once: it is what the process RUNS,
#: whatever the files on disk say later.
RUNNING: dict = {**_describe(), "started": time.strftime("%Y-%m-%d %H:%M"),
                 "dirty": bool(_git("status", "--porcelain"))}

_lock = threading.Lock()
_remote: dict = {"checked": 0.0, "commit": None, "error": None}


def _remote_main() -> dict:
    """origin/main's commit, cached for CHECK_SECONDS."""
    with _lock:
        if time.time() - _remote["checked"] < CHECK_SECONDS:
            return dict(_remote)
    out = _git("ls-remote", "origin", "refs/heads/main")
    with _lock:
        _remote["checked"] = time.time()
        if out:
            _remote["commit"], _remote["error"] = out.split()[0], None
        else:
            _remote["error"] = "GitHub could not be reached"
        return dict(_remote)


def status() -> dict:
    """Everything the badge shows.

    ``state`` is one of:
      "current"  running = on disk = GitHub's main
      "restart"  the checkout on disk is newer than what is running
      "behind"   GitHub's main is newer than the checkout
      "unknown"  not a git checkout, or GitHub not reachable
    """
    on_disk = _describe()
    remote = _remote_main()
    running = RUNNING.get("commit")
    if not running:
        state = "unknown"
    elif on_disk.get("commit") and on_disk["commit"] != running:
        state = "restart"
    elif remote.get("commit") and remote["commit"] != running:
        state = "behind"
    elif remote.get("error"):
        state = "unknown"
    else:
        state = "current"
    return {"state": state, "running": RUNNING, "on_disk": on_disk,
            "remote_commit": remote.get("commit"),
            "remote_error": remote.get("error")}
