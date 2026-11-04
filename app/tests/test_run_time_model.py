"""Size-based run time model. Run from the repo root:
python -m unittest app.tests.test_run_time_model"""

import shutil
import tempfile
import time
import unittest
from pathlib import Path

from app.services import phase_history, run_time_model, runner

J445 = {"cells": 2857074, "neighbors": 15, "threads": 27, "sgd_threads": 1}
TOTALS = {"load": 3, "extract": 354, "embed": 200, "plot": 34}


#: Complete runs logged on highmem2 (2026-11-03/04): size, and measured
#: minutes for the whole run. The model's constants were fitted on these and
#: the other logged runs, so this checks the formulas are wired correctly,
#: not how well they predict new data (see run_time_model's docstring).
LOGGED = [
    # fit, k, plots, supercell build skipped, measured minutes
    (98526, 30, 7, 12.3),
    (500000, 15, 7, 14.5),
    (500000, 15, 34, 19.6),
    (2857074, 15, 34, 92.6),
    (14876, 15, 34, 11.5),
]


class Predict(unittest.TestCase):
    def test_logged_j445_runs_within_20_percent(self):
        for fit, k, plots, minutes in LOGGED:
            info = dict(J445, fit_cells=fit, neighbors=k, spectral=1)
            totals = {"extract": 354, "embed": 200, "plot": plots}
            p = run_time_model.predict(info, totals)
            with self.subTest(fit=fit, plots=plots):
                self.assertAlmostEqual(sum(p.values()) / 60, minutes,
                                       delta=0.2 * minutes)

    def test_fit_on_all_cells_is_dominated_by_the_initialisation(self):
        p = run_time_model.predict(dict(J445, fit_cells=2857074, spectral=1),
                                   TOTALS)
        self.assertGreater(p["embed_init"], 0.5 * sum(p.values()))
        pca = run_time_model.predict(dict(J445, fit_cells=2857074, spectral=0),
                                     TOTALS)
        self.assertLess(pca["embed_init"], 120)

    def test_more_threads_less_time(self):
        few = run_time_model.predict(dict(J445, fit_cells=300000, threads=6),
                                     TOTALS)
        many = run_time_model.predict(dict(J445, fit_cells=300000), TOTALS)
        self.assertGreater(few["embed_nn"], 1.5 * many["embed_nn"])
        # The spectral initialisation runs on one thread.
        self.assertEqual(few["embed_init"], many["embed_init"])

    def test_nothing_before_the_size_is_known(self):
        self.assertIsNone(run_time_model.predict({}, {"load": 3}))


class Calibrate(unittest.TestCase):
    def test_host_twice_as_slow(self):
        info = dict(J445, fit_cells=300000)
        guess = run_time_model.predict(info, TOTALS)
        entry = {"phases": {"embed_init": 2 * guess["embed_init"]},
                 "info": info, "totals": TOTALS}
        fixed = run_time_model.calibrate(guess, [entry])
        self.assertAlmostEqual(fixed["embed_init"], 2 * guess["embed_init"])
        self.assertEqual(fixed["embed_nn"], guess["embed_nn"])


class FirstRunGetsAnEstimate(unittest.TestCase):
    """No run with these settings yet: the ::runinfo:: line gives a time."""

    def setUp(self):
        self.out = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.out, ignore_errors=True)
        cache = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, cache, ignore_errors=True)
        old = phase_history.HISTORY_FILE
        phase_history.HISTORY_FILE = Path(cache) / "h.json"
        self.addCleanup(setattr, phase_history, "HISTORY_FILE", old)

    def test_model_estimate_and_hidden_line(self):
        lines = [
            "::progress:: phase=load done=0 total=3 elapsed=1.0 label=a",
            "::progress:: phase=extract done=354 total=354 elapsed=2.0 label=b",
            "::runinfo:: cells=2857074 fit_cells=300000 neighbors=15 "
            "threads=27 sgd_threads=1",
            "::progress:: phase=embed_nn done=0 total=3 elapsed=3.0 label=c",
            "visible line",
        ]
        script = "; ".join(f"echo '{ln}'" for ln in lines) + "; sleep 3"
        r = runner.Runner()
        ok, _ = r.start(["/bin/sh", "-c", script], cwd="/", kind="run")
        self.assertTrue(ok)
        for _ in range(60):
            if "visible line" in r.snapshot()["log"]:
                break
            time.sleep(0.05)
        snap = r.snapshot()
        self.assertNotIn("::runinfo::", snap["log"])
        view = snap["progress"]
        self.assertEqual(view["estimate_kind"], "model")
        # 38 s neighbours + 31 s layout + 368 s optimisation and projection
        # + 149 s plots.
        self.assertAlmostEqual(view["remaining"], 38 + 31 + 368 + 149,
                               delta=60)
        r.cancel()


if __name__ == "__main__":
    unittest.main()
