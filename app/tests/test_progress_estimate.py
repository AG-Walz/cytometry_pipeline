"""Progress estimate. Run from the repo root:
python -m unittest app.tests.test_progress_estimate"""

import os
import shutil
import tempfile
import time
import unittest
from pathlib import Path

from app.services import detached, phase_history, progress_estimate, runner

#: Phase durations of the J445 v5 run of 2026-11-03 (93 min in total).
J445 = {"load": 106.0, "extract": 112.0, "embed_nn": 543.0,
        "embed_init": 3115.0, "embed": 1544.0, "plot": 171.0}


def _line(phase, done, total, elapsed, label="x"):
    return (f"::progress:: phase={phase} done={done} total={total} "
            f"elapsed={elapsed:.1f} label={label}")


def _complete_run(start=10.0):
    """Progress lines of one complete run with the J445 durations."""
    t, out = start, []
    for phase in progress_estimate.PHASE_ORDER:
        out.append(_line(phase, 0, 3 if phase != "embed_init" else 1, t))
        t += J445[phase]
        out.append(_line(phase, 3 if phase != "embed_init" else 1,
                         3 if phase != "embed_init" else 1, t))
    out.append("Done. Plots written to: /x")
    return out


def _until_init(start=10.0):
    """A run that has done load, extract and the neighbour search."""
    t, out = start, []
    for phase in ("load", "extract", "embed_nn"):
        out.append(_line(phase, 0, 3, t))
        t += J445[phase]
    out.append(_line("embed_nn", 3, 3, t))
    out.append(_line("embed_init", 0, 1, t))
    return out, t


class Durations(unittest.TestCase):
    def test_start_to_start(self):
        times = {"load": [10.0, 50.0], "extract": [80.0, 100.0],
                 "embed_nn": [120.0, 200.0]}
        self.assertEqual(progress_estimate.durations(times),
                         {"load": 70.0, "extract": 40.0, "embed_nn": 80.0})


class Estimate(unittest.TestCase):
    def _times(self):
        times = {}
        lines, t_init = _until_init()
        for ln in lines:
            m = runner._PROGRESS_RE.match(ln)
            runner._note_time(times, m.group(1), m.group(4))
        return times, t_init

    def test_one_step_phase_counts_down(self):
        """In a 0 / 1 phase the estimate goes DOWN as the clock runs."""
        times, t0 = self._times()
        progress = {"phase": "embed_init", "done": 0, "total": 1}
        left = [progress_estimate.estimate(progress, times, set(), J445,
                                           t0 + dt)["remaining"]
                for dt in (0, 600, 1800, 3000)]
        self.assertEqual(left, sorted(left, reverse=True))
        self.assertAlmostEqual(left[0], J445["embed_init"] + J445["embed"]
                               + J445["plot"], delta=60)
        self.assertAlmostEqual(left[1] - left[2], 1200, delta=1)

    def test_overrun(self):
        times, t0 = self._times()
        progress = {"phase": "embed_init", "done": 0, "total": 1}
        view = progress_estimate.estimate(progress, times, set(), J445,
                                          t0 + 4000)
        self.assertTrue(view["overrun"])
        self.assertAlmostEqual(view["remaining"], J445["embed"] + J445["plot"],
                               delta=1)

    def test_no_reference_no_number(self):
        times, t0 = self._times()
        progress = {"phase": "embed_init", "done": 0, "total": 1}
        view = progress_estimate.estimate(progress, times, set(), None,
                                          t0 + 60)
        self.assertIsNone(view["remaining"])
        self.assertFalse(view["has_reference"])
        self.assertEqual(len(view["phases"]), 6)

    def test_reused_phases_not_counted(self):
        times = {}
        for ln in (_line("load", 0, 3, 10), _line("extract", 0, 3, 40),
                   _line("embed_nn", 0, 3, 60), _line("plot", 1, 3, 61)):
            m = runner._PROGRESS_RE.match(ln)
            runner._note_time(times, m.group(1), m.group(4))
        resumed = progress_estimate.resumed_phases(
            ["Resuming: 'embedding' reused from checkpoints/embedding.rds"])
        view = progress_estimate.estimate(
            {"phase": "plot", "done": 1, "total": 3}, times, resumed, J445,
            62)
        states = {p["phase"]: p["state"] for p in view["phases"]}
        self.assertEqual(states["embed_init"], "reused")
        self.assertEqual(states["embed"], "reused")
        self.assertLess(view["remaining"], J445["plot"] * 1.6)


class History(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.dir, ignore_errors=True)
        self.path = Path(self.dir) / "h.json"

    def test_key_ignores_output_folder_and_resume(self):
        a = ["Rscript", "x.R", "--wsp", "w", "--output_dir", "/a"]
        b = ["Rscript", "x.R", "--wsp", "w", "--output_dir", "/b", "--no_resume"]
        c = ["Rscript", "x.R", "--wsp", "other", "--output_dir", "/a"]
        self.assertEqual(phase_history.key_for(a), phase_history.key_for(b))
        self.assertNotEqual(phase_history.key_for(a), phase_history.key_for(c))

    def test_save_merges(self):
        argv = ["Rscript", "x.R"]
        phase_history.save(argv, {"load": 10, "embed": 50}, self.path)
        phase_history.save(argv, {"load": 12}, self.path)
        self.assertEqual(phase_history.load(argv, self.path),
                         {"load": 12.0, "embed": 50.0})


class LogWithSeveralRuns(unittest.TestCase):
    """The 16 h 31 min bug: an output folder's log holds every run into it."""

    def setUp(self):
        self.out = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.out, ignore_errors=True)
        self.cache = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.cache, ignore_errors=True)
        old = phase_history.HISTORY_FILE
        phase_history.HISTORY_FILE = Path(self.cache) / "h.json"
        self.addCleanup(setattr, phase_history, "HISTORY_FILE", old)

    def test_only_last_run_counts(self):
        current, t_init = _until_init()
        text = "\n".join(
            [detached.RUN_MARKER + "2026-11-03 14:00:00 ──"] + _complete_run()
            + [detached.RUN_MARKER + "2026-11-04 10:15:00 ──"] + current) + "\n"
        log_path, status_path = detached.paths_for(self.out)
        log_path.parent.mkdir(parents=True)
        log_path.write_text(text, encoding="utf-8")
        now = time.time()
        status_path.write_text(
            '{"pid": %d, "argv": ["Rscript", "x.R"], "kind": "run", '
            '"started_at": %f, "running": true}' % (os.getpid(), now - t_init - 600))
        real = detached.is_alive
        detached.is_alive = lambda pid: True
        self.addCleanup(setattr, detached, "is_alive", real)

        r = runner.Runner()
        self.assertTrue(r.attach(self.out))
        view = r.snapshot()["progress"]
        self.assertEqual(view["phase"], "embed_init")
        self.assertTrue(view["has_reference"])
        expected = (J445["embed_init"] - 600) + J445["embed"] + J445["plot"]
        self.assertAlmostEqual(view["remaining"], expected,
                               delta=0.1 * expected)
        states = {p["phase"]: p["state"] for p in view["phases"]}
        self.assertEqual(states["embed"], "pending")
        self.assertEqual(states["plot"], "pending")


class RecordsHistory(unittest.TestCase):
    def setUp(self):
        self.out = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.out, ignore_errors=True)
        self.cache = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.cache, ignore_errors=True)
        old = phase_history.HISTORY_FILE
        phase_history.HISTORY_FILE = Path(self.cache) / "h.json"
        self.addCleanup(setattr, phase_history, "HISTORY_FILE", old)

    def test_complete_detached_run_is_saved(self):
        script = "; ".join(f"echo '{ln}'" for ln in (
            _line("load", 0, 1, 1), _line("extract", 0, 2, 5),
            _line("extract", 2, 2, 9), "Done. Plots written to: /x"))
        argv = ["/bin/sh", "-c", script]
        r = runner.Runner()
        ok, _ = r.start_detached(argv, cwd="/", output_dir=self.out, kind="run")
        self.assertTrue(ok)
        for _ in range(100):
            if r.snapshot()["status"] != runner.RUNNING:
                break
            time.sleep(0.05)
        self.assertEqual(r.snapshot()["status"], runner.DONE)
        self.assertEqual(phase_history.load(argv),
                         {"load": 4.0, "extract": 4.0})


class Bar(unittest.TestCase):
    """The bar is one Mantine progress bar; the time is for the whole run."""

    def test_sections_and_whole_run_time(self):
        from app.callbacks import execution
        times, t0 = Estimate._times(self)
        view = progress_estimate.estimate(
            {"phase": "embed_init", "done": 0, "total": 1, "label": "x"},
            times, set(), J445, t0 + 600)
        bar = execution._phase_bar(view)
        self.assertEqual(type(bar).__name__, "ProgressRoot")
        self.assertAlmostEqual(sum(sec.value for sec in bar.children), 100,
                               delta=1)
        running = [sec for sec in bar.children if getattr(sec, "striped", False)]
        self.assertEqual(len(running), 1)
        # 42 min of this step + 26 min of optimising + 3 min of plots.
        self.assertEqual(execution._format_remaining(view),
                         "about 1 h 10 min left")


if __name__ == "__main__":
    unittest.main()
