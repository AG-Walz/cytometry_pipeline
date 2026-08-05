"""Runs that outlive the app.

A UMAP over 7.9 million cells is an overnight job. The app was never in a
position to host one: the process was a child of gunicorn with its output going
into a pipe and its state in a module-level object, so restarting the app — or
deploying to it, or having it crash — left a running R process nobody could see
and a log nobody could read.

Nothing here supervises the process. The child is already immune to the app
dying: ``start_new_session=True`` puts it in its own session, and losing a
parent does not kill a session leader. What was missing was the ability to FIND
it again, so this writes down the two things that makes possible:

    <output_dir>/_run/run.log      everything the process printed
    <output_dir>/_run/status.json  pid, argv, when it started, what kind

The log goes to a FILE rather than a pipe for the same reason. A pipe only
exists while both ends do; a file is still there in the morning, and a second
reader — a person with `tail -f`, a second browser — costs nothing.

Reattaching is then: read status.json, check whether that pid is still alive,
and tail the log from the beginning. That works for a run this app started, a
run the PREVIOUS app instance started, and a run someone launched by hand with
the command the UI printed — which is the point of printing it.
"""

from __future__ import annotations

import json
import os
import signal
import subprocess
import time
from pathlib import Path

#: Subdirectory of the run's output directory holding its bookkeeping. Inside
#: the output directory on purpose: the run's log belongs with the run's
#: results, and deleting the results should take the log with them.
RUN_SUBDIR = "_run"

LOG_NAME = "run.log"
STATUS_NAME = "status.json"


def paths_for(output_dir: str | os.PathLike) -> tuple[Path, Path]:
    d = Path(output_dir) / RUN_SUBDIR
    return d / LOG_NAME, d / STATUS_NAME


def is_alive(pid: int) -> bool:
    """Whether ``pid`` is still DOING something.

    signal 0 checks for existence without delivering anything, and a pid we are
    not permitted to signal still exists — a run started by another user is
    running whatever we are allowed to do to it.

    But existence is not the question. A run this app started is its child, and
    nothing ever wait()s for it: when it exits it becomes a ZOMBIE, which
    exists, answers signal 0, and will keep doing so until the app itself dies.
    Asking only "does this pid exist" therefore reported every finished run as
    still running, forever. The process state settles it — and the check is
    skipped entirely on a run started by a previous app instance, which was
    re-parented to init and reaped properly, so this costs nothing in the case
    that already worked.
    """
    if not pid or pid <= 0:
        return False
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    except OSError:
        return False

    # Field 3 of /proc/<pid>/stat is the state code. The command name in field
    # 2 is wrapped in parentheses and may itself contain spaces and brackets,
    # so the split is anchored on the LAST ')' rather than on whitespace.
    try:
        with open(f"/proc/{pid}/stat", "r", encoding="utf-8") as fh:
            after_name = fh.read().rsplit(")", 1)[1].split()
        if after_name and after_name[0] == "Z":
            return False
    except (OSError, IndexError):
        pass                       # not Linux, or it vanished — signal 0 stands
    return True


def launch(argv: list[str], cwd: str, output_dir: str, kind: str,
           meta: dict | None = None) -> dict:
    """Start ``argv`` detached, writing its output into ``output_dir``.

    Returns the status dict that was written. Raises OSError if the process
    cannot be started or the directory cannot be created — the caller reports
    that; there is nothing sensible to do with a half-started run.
    """
    log_path, status_path = paths_for(output_dir)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    # Append, not truncate: a re-run into the same output directory should not
    # erase the evidence of the one before it, which is often exactly what the
    # person is trying to compare against.
    log_file = open(log_path, "a", buffering=1, encoding="utf-8",
                    errors="replace")
    try:
        log_file.write(f"\n── started {time.strftime('%Y-%m-%d %H:%M:%S')} ──\n")
        proc = subprocess.Popen(
            argv,
            cwd=cwd,
            stdout=log_file,
            stderr=subprocess.STDOUT,
            stdin=subprocess.DEVNULL,
            text=True,
            # Its own session: the app can exit, be redeployed or crash without
            # taking the run with it. This is the whole feature.
            start_new_session=True,
        )
    finally:
        # The child holds its own duplicate of the descriptor; this one is not
        # needed and holding it would keep the file open for the app's lifetime.
        log_file.close()

    status = {
        "pid": proc.pid,
        "argv": list(argv),
        "cwd": cwd,
        "kind": kind,
        "output_dir": str(output_dir),
        "started_at": time.time(),
        "meta": dict(meta or {}),
    }
    status_path.write_text(json.dumps(status, indent=1), encoding="utf-8")
    return status


def read_status(output_dir: str | os.PathLike) -> dict | None:
    """The recorded status of a run in ``output_dir``, or None."""
    _, status_path = paths_for(output_dir)
    try:
        status = json.loads(status_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    if not isinstance(status, dict) or "pid" not in status:
        return None
    status["running"] = is_alive(int(status.get("pid") or 0))
    return status


def read_log(output_dir: str | os.PathLike, max_bytes: int = 4_000_000) -> str:
    """The run's log, tail-truncated.

    An overnight run over 2,187 plots writes a lot; the whole file is not worth
    shipping into a browser on every poll. The tail is what anyone wants — and
    the truncation says so rather than silently starting mid-sentence.
    """
    log_path, _ = paths_for(output_dir)
    try:
        size = log_path.stat().st_size
        with log_path.open("r", encoding="utf-8", errors="replace") as fh:
            if size > max_bytes:
                fh.seek(size - max_bytes)
                fh.readline()          # drop the partial line we landed inside
                head = (f"… earlier output omitted "
                        f"({(size - max_bytes) // 1024} KB); the full log is at "
                        f"{log_path}\n")
                return head + fh.read()
            return fh.read()
    except OSError:
        return ""


def stop(output_dir: str | os.PathLike) -> bool:
    """Terminate a detached run. Returns whether a signal was delivered.

    The whole process GROUP, because the run forks plotting workers and killing
    only the leader would leave them running with nothing to report to.
    """
    status = read_status(output_dir)
    if not status or not status.get("running"):
        return False
    pid = int(status["pid"])
    try:
        os.killpg(os.getpgid(pid), signal.SIGTERM)
        return True
    except (ProcessLookupError, PermissionError, OSError):
        try:
            os.kill(pid, signal.SIGTERM)
            return True
        except OSError:
            return False


def find_recent(parents: list[str], limit: int = 20) -> list[dict]:
    """Runs found under the given parent directories, newest first.

    Used to offer "you have a run going" after the app restarts. Only one level
    deep: the app writes each run into its own subfolder of the output parent,
    so a recursive walk would buy nothing and could wander into a data tree.
    """
    found: list[dict] = []
    for parent in parents:
        p = Path(parent or "")
        if not p.is_dir():
            continue
        try:
            entries = [d for d in p.iterdir() if d.is_dir()]
        except OSError:
            continue
        for d in entries:
            status = read_status(d)
            if status:
                found.append(status)
    found.sort(key=lambda s: s.get("started_at") or 0, reverse=True)
    return found[:limit]
