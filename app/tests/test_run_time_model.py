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


class Predict(unittest.TestCase):
    def test_reproduces_the_j445_fit_on_all_cells(self):
        """Measured on highmem2: 543 s, 3115 s, 1544 s (93 min in total)."""
        p = run_time_model.predict(dict(J445, fit_cells=2857074), TOTALS)
        self.assertAlmostEqual(p["embed_nn"], 543, delta=30)
        self.assertAlmostEqual(p["embed_init"], 3115, delta=150)
        self.assertAlmostEqual(p["embed"], 1544, delta=80)
        self.assertAlmostEqual(sum(p.values()) / 60, 93, delta=6)

    def test_reproduces_the_j445_fit_on_300k_cells(self):
        """Measured: 73 s, 164 s, 396 s (optimisation and projection)."""
        p = run_time_model.predict(dict(J445, fit_cells=300000), TOTALS)
        self.assertAlmostEqual(p["embed_nn"], 73, delta=15)
        self.assertAlmostEqual(p["embed_init"], 164, delta=30)
        self.assertAlmostEqual(p["embed"], 396, delta=40)

    def test_more_threads_less_time(self):
        few = run_time_model.predict(dict(J445, fit_cells=300000, threads=6),
                                     TOTALS)
        many = run_time_model.predict(dict(J445, fit_cells=300000), TOTALS)
        self.assertGreater(few["embed_nn"], 2 * many["embed_nn"])
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
        # 71 s neighbours + 166 s layout + 396 s optimisation + plots.
        self.assertAlmostEqual(view["remaining"], 71 + 166 + 396 + 171,
                               delta=60)
        r.cancel()


if __name__ == "__main__":
    unittest.main()
