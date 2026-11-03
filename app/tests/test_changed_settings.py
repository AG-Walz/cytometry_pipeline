"""Changed-from-defaults summary. Run from the repo root:
python -m unittest app.tests.test_changed_settings"""

import unittest

from app import config
from app.callbacks import common


def form(**over):
    f = {**config.DEFAULTS, **config.FORM_START}
    f.update(over)
    return f


class Changed(unittest.TestCase):
    def labels(self, **over):
        return [label for label, _, _ in common.changed_settings(form(**over))]

    def test_defaults_report_nothing(self):
        self.assertEqual(self.labels(), [])

    def test_changed_values_are_named(self):
        self.assertEqual(self.labels(n_neighbors=80, point_size=0.1),
                         ["n_neighbors", "Point size"])

    def test_number_types_do_not_count_as_changes(self):
        self.assertEqual(self.labels(n_neighbors=15.0, output_formats=("png",)), [])

    def test_min_dist_default_follows_the_layout(self):
        self.assertEqual(self.labels(layout="gates", min_dist=0.25), ["Layout"])

    def test_gate_settings_only_count_in_the_gate_layout(self):
        self.assertNotIn("Gap at each gate cut", self.labels(gate_gap=3))
        self.assertIn("Gap at each gate cut", self.labels(layout="gates", gate_gap=3))


class Preview(unittest.TestCase):
    def test_preview_is_small_and_separate(self):
        from app.services import rscript
        params = common.run_params(form(wsp="/w.wsp", fcs_dir="/f",
                                        output_parent="/o", output_name="run1",
                                        n_neighbors=30), ["/P"], {})
        p = common.preview_params(params)
        self.assertTrue(p["output_dir"].endswith("_preview"))
        a = rscript.run_argv(p)
        self.assertEqual(a[a.index("--downsample") + 1], "1000")
        self.assertEqual(a[a.index("--n_neighbors") + 1], "30")   # settings kept
        self.assertNotEqual(params["output_dir"], p["output_dir"])


class LoadSettings(unittest.TestCase):
    def test_a_recorded_run_becomes_form_values(self):
        opts = {"layout": "markers", "normalize_data": "zscore",
                "n_neighbors": 30, "min_dist": 0.05, "point_size": None,
                "output_formats": "png,svg", "no_resume": True,
                "umap_subsample": -1, "clip_quantile": 0.001,
                "gate_separation": 0.75, "wsp": "/somewhere.wsp"}
        f = common.form_from_settings(opts)
        self.assertEqual(f["normalize"], "zscore")
        self.assertEqual(f["n_neighbors"], 30)
        self.assertIsNone(f["point_size"])
        self.assertEqual(f["output_formats"], ["png", "svg"])
        self.assertFalse(f["resume"])
        self.assertEqual(f["random_fit_cells"], config.DEFAULTS["random_fit_cells"])
        self.assertEqual(f["clip_quantile"], "0.001")
        self.assertNotIn("wsp", f)        # inputs are not settings

    def test_gates_fit_budget(self):
        f = common.form_from_settings({"layout": "gates", "umap_subsample": -1})
        self.assertEqual(f["umap_subsample"], 0)


if __name__ == "__main__":
    unittest.main()
