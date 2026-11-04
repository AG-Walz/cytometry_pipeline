"""A single-slot subprocess runner with a live, pollable log.

The app is deliberately one-run-at-a-time, so one module-level ``Runner``
instance owns whichever R process is active — an inspect or a full run. Both
stream into the same log, and neither can start while the other is alive, which
is exactly the guarantee the UI promises.

Output is read on a background thread into a bounded buffer; the UI polls
``snapshot()`` from a ``dcc.Interval``. Nothing here is safe for multi-user
serving, and it is not meant to be: the app binds to localhost.
"""

from __future__ import annotations

import os
import re
import signal
import subprocess
import threading
import time

from app.services import detached, phase_history, progress_estimate, run_time_model

#: Hard cap on retained log lines, so a pathological run cannot exhaust memory.
MAX_LINES = 20000

#: Progress lines emitted by ``emit_progress()`` in src/R_shared/io_utils.R.
#: They are intercepted here rather than shown: they are for the progress bar,
#: and one per sample would bury the log they share a stream with. R owns the
#: counting; this only reads it. Keep in step with PROGRESS_PREFIX there.
#: ``elapsed`` is optional so a log written before R started emitting it still
#: parses — the bar then loses its timings, not its position. It sits before
#: ``label``, which stays last because a label may contain anything at all.
_PROGRESS_RE = re.compile(
    r"^::progress:: phase=(\S+) done=(\d+) total=(\d+)"
    r"(?: elapsed=([0-9.]+))? label=(.*)$"
)

#: The size of the job, from emit_run_info() in src/R_shared/io_utils.R:
#: ``::runinfo:: cells=2857074 fit_cells=300000 neighbors=15 threads=27``.
#: Hidden from the log, like the progress lines.
_RUNINFO_RE = re.compile(r"^::runinfo:: (.*)$")


def _parse_run_info(text: str) -> dict[str, float]:
    out = {}
    for pair in text.split():
        key, _, value = pair.partition("=")
        try:
            out[key] = float(value)
        except ValueError:
            continue
    return out


#: The phases, their labels and order live with the estimate that uses them.
#: Re-exported here because callers have always found them on the runner.
PROGRESS_PHASES = progress_estimate.PROGRESS_PHASES
PHASE_WEIGHTS = progress_estimate.PHASE_WEIGHTS
PHASE_ORDER = progress_estimate.PHASE_ORDER

IDLE = "idle"
RUNNING = "running"
DONE = "done"
FAILED = "failed"
CANCELLED = "cancelled"


def _note_time(times: dict[str, list[float]], phase: str,
               elapsed: str | None) -> None:
    """Record when a phase was first and last heard from.

    ``elapsed`` is R's own seconds-since-start (see emit_progress in
    src/R_shared/io_utils.R), so it needs no shared clock and survives the app
    being restarted mid-run. A line without it — an older R, or another
    emitter — is simply not counted.
    """
    if elapsed is None:
        return
    try:
        value = float(elapsed)
    except ValueError:
        return
    span = times.get(phase)
    if span is None:
        times[phase] = [value, value]
    else:
        span[1] = value


def _reference_from_runs(runs: list[dict]) -> dict[str, float] | None:
    """Phase durations of the earlier COMPLETE runs in one log, newest first.

    A run counts as complete when it printed the line R prints last on
    success. Reused phases took no time, so an older run fills them in.
    """
    out: dict[str, float] = {}
    for run in reversed(runs):
        if not run["times"] or not Runner._looks_finished(run["lines"]):
            continue
        resumed = progress_estimate.resumed_phases(run["lines"])
        for phase, took in progress_estimate.durations(run["times"]).items():
            if phase not in resumed and phase not in out:
                out[phase] = took
    return out or None


class Runner:
    def __init__(self) -> None:
        self._lock = threading.Lock()
        self._proc: subprocess.Popen | None = None
        self._thread: threading.Thread | None = None
        self._lines: list[str] = []
        self._truncated = False
        self._status = IDLE
        self._kind = ""
        self._argv: list[str] = []
        self._meta: dict = {}
        self._returncode: int | None = None
        self._started_at: float | None = None
        self._finished_at: float | None = None
        self._cancelling = False
        #: Bumped on every completion so the UI can act once per finished run.
        self._generation = 0
        #: Latest ::progress:: line of the CURRENT run.
        self._progress: dict = {}
        #: {phase: [first_elapsed, last_elapsed]} in the R process's own
        #: seconds, for the current run only.
        self._phase_times: dict[str, list[float]] = {}
        #: Phases the current run skips because it reused a saved step.
        self._resumed: set[str] = set()
        #: Phase durations of the reference run, or None; see
        #: app/services/progress_estimate.py.
        self._reference: dict[str, float] | None = None
        #: The size of the current run (its ::runinfo:: line), the latest
        #: total of each phase, and the saved runs the model calibrates on.
        self._run_info: dict[str, float] = {}
        self._phase_totals: dict[str, int] = {}
        self._history_entries: list[dict] = []
        #: The detached run whose durations are already in the history.
        self._recorded_dir: str | None = None
        #: Set when following a run that is not our child (see attach()).
        self._detached_dir: str | None = None
        self._detached_pid: int = 0

    # -- lifecycle ---------------------------------------------------------

    def start(self, argv: list[str], cwd: str, kind: str, meta: dict | None = None):
        """Launch ``argv``. Returns ``(ok, message)``; refuses while busy."""
        with self._lock:
            if self._status == RUNNING:
                return False, "A run is already in progress."
            self._lines = []
            self._truncated = False
            self._status = RUNNING
            self._kind = kind
            self._argv = list(argv)
            self._meta = dict(meta or {})
            self._returncode = None
            self._progress = {}
            self._phase_times = {}
            self._resumed = set()
            self._reference = phase_history.load(argv)
            self._run_info, self._phase_totals = {}, {}
            self._history_entries = phase_history.entries()
            self._started_at = time.time()
            self._finished_at = None
            self._cancelling = False
            # This run's output comes down a pipe. Any detached run we were
            # following is no longer what the UI is showing, and leaving the
            # pointer set made snapshot() re-read THAT run's log on every poll —
            # overwriting this one's lines and re-deciding this one's status from
            # a finished process's pid. That is what surfaced as an inspect
            # reporting "failed (exit code None)" while printing a complete log.
            self._detached_dir = None
            self._detached_pid = 0

            try:
                # start_new_session puts the child in its own process group so
                # cancelling kills R and anything it spawned, not just the
                # stdbuf wrapper.
                self._proc = subprocess.Popen(
                    argv,
                    cwd=cwd,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                    bufsize=1,
                    start_new_session=True,
                )
            except (OSError, ValueError) as exc:
                self._status = FAILED
                self._finished_at = time.time()
                self._generation += 1
                self._lines = [f"Failed to launch: {exc}", "", "Command was:"]
                self._lines += ["  " + " ".join(argv)]
                return False, str(exc)

            self._thread = threading.Thread(
                target=self._pump, args=(self._proc,), daemon=True
            )
            self._thread.start()
            return True, "started"

    def _pump(self, proc: subprocess.Popen) -> None:
        """Drain the child's output, then record its exit status."""
        assert proc.stdout is not None
        try:
            for line in proc.stdout:
                line = line.rstrip("\n")
                hit = _PROGRESS_RE.match(line)
                if hit:
                    phase, done, total, elapsed, label = hit.groups()
                    with self._lock:
                        self._progress = {
                            "phase": phase,
                            "done": int(done),
                            "total": int(total),
                            "label": label,
                        }
                        self._note_phase_time(phase, elapsed)
                        self._phase_totals[phase] = int(total)
                    continue
                info = _RUNINFO_RE.match(line)
                if info:
                    with self._lock:
                        self._run_info = _parse_run_info(info.group(1))
                    continue
                with self._lock:
                    self._resumed |= progress_estimate.resumed_phases([line])
                    if len(self._lines) < MAX_LINES:
                        self._lines.append(line)
                    elif not self._truncated:
                        self._truncated = True
                        self._lines.append(
                            f"… output truncated at {MAX_LINES} lines; "
                            "the run continues and files are still being written."
                        )
        except (ValueError, OSError):
            # Pipe closed underneath us (e.g. during a cancel).
            pass
        finally:
            returncode = proc.wait()
            with self._lock:
                self._returncode = returncode
                self._finished_at = time.time()
                if self._cancelling:
                    self._status = CANCELLED
                elif returncode == 0:
                    self._status = DONE
                else:
                    self._status = FAILED
                self._record_history()
                self._generation += 1
                self._proc = None

    def cancel(self) -> bool:
        """Terminate the active process group. No-op when idle."""
        if self._detached_dir:
            stopped = detached.stop(self._detached_dir)
            if stopped:
                with self._lock:
                    self._cancelling = True
                    self._status = CANCELLED
                    self._finished_at = time.time()
                    self._generation += 1
            return stopped
        with self._lock:
            proc = self._proc
            if proc is None or self._status != RUNNING:
                return False
            self._cancelling = True
            self._lines.append("")
            self._lines.append("── cancelled by user ──")
            pid = proc.pid

        try:
            os.killpg(os.getpgid(pid), signal.SIGTERM)
        except (ProcessLookupError, PermissionError, OSError):
            try:
                proc.terminate()
            except OSError:
                return False

        def _hard_kill() -> None:
            time.sleep(5)
            with self._lock:
                still = self._proc
            if still is not None and still.pid == pid:
                try:
                    os.killpg(os.getpgid(pid), signal.SIGKILL)
                except OSError:
                    pass

        threading.Thread(target=_hard_kill, daemon=True).start()
        return True

    # -- detached runs -----------------------------------------------------

    def start_detached(self, argv, cwd, output_dir, kind, meta=None):
        """Launch a run that outlives the app, and follow its log.

        Same single-slot rule as start(): one thing at a time. The difference is
        where the output goes — a file in the run's own output directory rather
        than a pipe — which is what makes it survivable and re-readable. See
        app/services/detached.py.
        """
        with self._lock:
            if self._status == RUNNING:
                return False, "A run is already in progress."
        try:
            status = detached.launch(argv, cwd, output_dir, kind, meta)
        except OSError as exc:
            with self._lock:
                self._status = FAILED
                self._lines = [f"Failed to launch: {exc}", "", "Command was:",
                               "  " + " ".join(argv)]
                self._finished_at = time.time()
                self._generation += 1
            return False, str(exc)
        self.attach(output_dir)
        return True, "started"

    def attach(self, output_dir) -> bool:
        """Follow a detached run — ours, or one a previous app instance began.

        Returns False when there is nothing there to follow. Refuses while
        something else is running, so attaching cannot silently displace the run
        the user is currently watching.
        """
        status = detached.read_status(output_dir)
        if not status:
            return False
        with self._lock:
            if self._status == RUNNING and self._detached_dir != str(output_dir):
                return False
            self._detached_dir = str(output_dir)
            self._argv = list(status.get("argv") or [])
            self._kind = status.get("kind") or "run"
            self._meta = dict(status.get("meta") or {})
            self._started_at = status.get("started_at")
            self._detached_pid = int(status.get("pid") or 0)
            self._returncode = status.get("exit_code")
            self._status = RUNNING if status.get("running") else DONE
            self._finished_at = None if status.get("running") else time.time()
            self._progress = {}
            self._phase_times = {}
            self._resumed = set()
            self._reference = phase_history.load(self._argv)
            self._run_info, self._phase_totals = {}, {}
            self._history_entries = phase_history.entries()
            self._recorded_dir = None
            self._lines = []
        self._refresh_detached()
        return True

    def _refresh_detached(self) -> None:
        """Re-read the detached run's log and re-derive status and progress.

        The log is the only channel — there is no pipe to a process we may not
        even have started — so progress is recovered by scanning it. Cheap
        enough at poll rate: the scan is a regex over the tail, and only the
        LAST progress line of each phase matters.
        """
        directory = self._detached_dir
        if not directory:
            return
        # Aliveness FIRST, then the log. The other order loses a race that
        # happens on every fast run: read a log that is still being written,
        # then observe a process that has since exited, and conclude from an
        # incomplete log that it stopped early. Asking about the process first
        # means a log read after "not alive" is a log nobody can still append
        # to.
        alive = detached.is_alive(self._detached_pid)
        text = detached.read_log(directory)
        code = None if alive else detached.read_exit_code(directory)
        lines = text.splitlines()

        # The log holds every run made into this output folder, one after the
        # other. Only the LAST run is the one in progress. Reading all of them
        # as one run is what made a 1 h 30 run show "16 h 31 min left".
        keep: list[str] = []
        runs: list[dict] = []
        for line in lines:
            if line.startswith(detached.RUN_MARKER) or not runs:
                runs.append({"progress": {}, "times": {}, "lines": [],
                             "info": {}, "totals": {}})
            hit = _PROGRESS_RE.match(line)
            if hit:
                phase, done, total, elapsed, label = hit.groups()
                runs[-1]["progress"] = {"phase": phase, "done": int(done),
                                        "total": int(total), "label": label}
                _note_time(runs[-1]["times"], phase, elapsed)
                runs[-1]["totals"][phase] = int(total)
                continue
            info = _RUNINFO_RE.match(line)
            if info:
                runs[-1]["info"] = _parse_run_info(info.group(1))
                continue
            runs[-1]["lines"].append(line)
            keep.append(line)
        current = runs[-1] if runs else {"progress": {}, "times": {},
                                          "lines": [], "info": {},
                                          "totals": {}}
        progress = current["progress"]
        resumed = progress_estimate.resumed_phases(current["lines"])
        log_reference = _reference_from_runs(runs[:-1])

        with self._lock:
            # The reads above ran without the lock, and start() may have begun
            # a NEW run in the meantime — a poll that started while a finished
            # run was still attached, then an inspect clicked a moment later.
            # What was read describes the old run; applying it here overwrote
            # the inspect's empty log with the old one and, since the old
            # process was not alive, declared the inspect finished (or failed
            # with no exit code and no log) a second after it started.
            if self._detached_dir != directory:
                return
            if len(keep) > MAX_LINES:
                self._truncated = True
                keep = keep[-MAX_LINES:]
            self._lines = keep
            # Rebuilt from the log every poll, because for a detached run the
            # log is the only memory there is: a browser reopened an hour into
            # a run has to recover how long each earlier phase took from what
            # R wrote, not from anything held here.
            self._phase_times = current["times"]
            self._resumed = resumed
            self._run_info = current["info"]
            self._phase_totals = current["totals"]
            self._progress = progress
            # A complete run with the same settings is the best reference. An
            # earlier complete run in this log is the next best: before this
            # history existed, it is the only one there is.
            history = phase_history.load(self._argv) if self._argv else None
            self._reference = history or log_reference
            if alive:
                return
            finished = code == 0 if code is not None else self._looks_finished(keep)
            verdict = DONE if finished else FAILED
            if self._status == RUNNING:
                self._returncode = code
                self._status = verdict
                self._finished_at = time.time()
                self._generation += 1
            elif self._status == FAILED and verdict == DONE:
                # A verdict taken from a log that turned out to be incomplete.
                # Correcting it matters more than never changing our mind: the
                # alternative is an app that says a finished run failed and has
                # no way to stop saying it.
                self._returncode = code
                self._status = DONE
                self._generation += 1
            if verdict == DONE and self._status == DONE:
                self._record_history_once(directory)

    def _note_phase_time(self, phase: str, elapsed: str | None) -> None:
        """Caller holds the lock. See _note_time."""
        _note_time(self._phase_times, phase, elapsed)

    def _record_history_once(self, directory: str) -> None:
        """Caller holds the lock. Once per run, however often it is polled.

        Also for a run that ended before the app attached to it (the app was
        restarted, or the run was very short), which never goes through the
        RUNNING -> DONE change.
        """
        if self._recorded_dir != directory:
            self._recorded_dir = directory
            self._record_history()

    @staticmethod
    def _looks_finished(lines: list[str]) -> bool:
        """Did the run reach its end, or did it stop?

        Only for runs started before the launcher recorded an exit code — see
        detached.read_exit_code(), which is the real answer when it exists.
        Completion is judged from what the run printed; the R side prints one
        of these last and only on success.
        """
        tail = "\n".join(lines[-40:])
        return "Done. Plots written to:" in tail or "Load report:" in tail

    # -- inspection --------------------------------------------------------

    def snapshot(self) -> dict:
        """A consistent view of the run for the UI."""
        # A detached run has no pipe feeding the buffer, so its log and status
        # are re-read here, at poll rate, before the view is taken.
        if self._detached_dir:
            self._refresh_detached()
        with self._lock:
            elapsed = None
            if self._started_at is not None:
                end = self._finished_at or time.time()
                elapsed = end - self._started_at
            return {
                "status": self._status,
                "kind": self._kind,
                "log": "\n".join(self._lines),
                "n_lines": len(self._lines),
                "returncode": self._returncode,
                "elapsed": elapsed,
                "argv": list(self._argv),
                "meta": dict(self._meta),
                "generation": self._generation,
                "progress": self._progress_view(),
            }

    def _progress_view(self) -> dict:
        """The bar and the time left; see app/services/progress_estimate.py.

        Caller holds the lock. Returns {} when nothing has reported yet, which
        is the honest state for the first minutes of a load: better a bar that
        appears late than one that sits at 0% claiming to know something.
        """
        p = self._progress
        if not p or not p.get("total"):
            return {}
        # The current time in R's seconds. R's own time stamps only move when
        # R prints a progress line, and a 50-minute step prints none, so the
        # clock comes from when the run started.
        run_now = (time.time() - self._started_at
                   if self._started_at and self._status == RUNNING else None)
        reference, kind = self._best_reference()
        view = progress_estimate.estimate(
            p, self._phase_times, self._resumed, reference, run_now, kind=kind)
        if self._status != RUNNING:
            view["remaining"] = None
        return view

    def _record_history(self) -> None:
        """Keep the phase durations of a run that finished, for the next one.

        Caller holds the lock. Reused phases took no time and are left out,
        and so is a run with no time stamps (an older R).
        """
        if self._status != DONE or not self._argv or not self._phase_times:
            return
        took = progress_estimate.durations(self._phase_times)
        took = {ph: v for ph, v in took.items() if ph not in self._resumed}
        phase_history.save(self._argv, took, info=self._run_info,
                           totals=self._phase_totals)

    def _best_reference(self) -> tuple[dict[str, float] | None, str]:
        """The reference for the estimate, and what kind it is.

        Caller holds the lock. A complete run with the same settings first;
        else a prediction from the size of the job, calibrated on this host's
        saved runs; else nothing (before R has said the size).
        """
        if self._reference:
            return self._reference, "measured"
        guess = run_time_model.predict(self._run_info, self._phase_totals)
        if guess:
            return run_time_model.calibrate(guess, self._history_entries), "model"
        return None, "measured"

    @property
    def busy(self) -> bool:
        with self._lock:
            return self._status == RUNNING


#: The app's single runner.
RUNNER = Runner()


def format_elapsed(seconds: float | None) -> str:
    if seconds is None:
        return "—"
    seconds = int(seconds)
    return f"{seconds // 60:02d}:{seconds % 60:02d}"
