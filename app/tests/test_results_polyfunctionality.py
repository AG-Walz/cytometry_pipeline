"""Polyfunctionality maps in the results gallery. Run from the repo root:
python -m unittest app.tests.test_results_polyfunctionality"""

import json
import tempfile
import unittest
from pathlib import Path

from app.services import results_scan

WHOLE = "umap_T_cells_polyfunctionality"
LEVELS = ("umap_T_cells_polyfunctionality_visit_V2", "umap_T_cells_polyfunctionality_visit_V0")


class Gallery(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.run = Path(self.tmp.name)
        (self.run / "overview").mkdir()
        (self.run / "overview" / "umap_T_cells_by_combination.png").write_bytes(b"x")
        (self.run / "overview" / "umap_T_cells_by_file.png").write_bytes(b"x")
        poly = self.run / "polyfunctionality"
        (poly / "visit").mkdir(parents=True)
        (poly / f"{WHOLE}.png").write_bytes(b"x")
        for stem in LEVELS:
            (poly / "visit" / f"{stem}.png").write_bytes(b"x")

    def test_indexed(self):
        idx = [{"file": f"polyfunctionality/{WHOLE}", "kind": "polyfunctionality",
                "column": None, "level": None, "reference": None}]
        idx += [{"file": f"polyfunctionality/visit/{s}", "kind": "polyfunctionality_split",
                 "column": "visit", "level": s.rsplit("_", 1)[1], "reference": None}
                for s in LEVELS]
        (self.run / "plot_index.json").write_text(json.dumps(idx))
        r = results_scan.scan(self.run)
        # Straight after the combination plot, before by_file.
        self.assertEqual([p.key for p in r.overview],
                         ["umap_T_cells_by_combination", WHOLE, "umap_T_cells_by_file"])
        # Per level under the split column, in the order R drew them.
        self.assertEqual([p.sample for p in r.split_poly["visit"]], ["V2", "V0"])
        self.assertIn("visit", r.split_columns)
        self.assertEqual(r.other, [])

    def test_without_index_by_name(self):
        r = results_scan.scan(self.run)
        self.assertIn(WHOLE, [p.key for p in r.overview])
        # Every level map lands among the polyfunctionality maps, none in "other".
        self.assertEqual(sorted(p.key for ps in r.split_poly.values() for p in ps),
                         sorted(LEVELS))
        self.assertEqual(r.other, [])
