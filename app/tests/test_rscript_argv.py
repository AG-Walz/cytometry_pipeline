"""The R command the app builds. Run from the repo root:
python -m unittest app.tests.test_rscript_argv"""

import unittest

from app import config
from app.services import rscript


def argv(**over):
    params = dict(config.DEFAULTS, wsp="/x.wsp", fcs_dir="/f", output_dir="/o",
                  population="/P", layout="markers")
    params.update(over)
    return rscript.run_argv(params)


class Normalisation(unittest.TestCase):
    def test_zscore_is_sent_as_normalize_data(self):
        a = argv(normalize="zscore")
        self.assertEqual(a[a.index("--normalize_data") + 1], "zscore")
        self.assertNotIn("--normalize", a)

    def test_default_sends_nothing(self):
        a = argv()
        self.assertNotIn("--normalize_data", a)
        self.assertNotIn("--normalize", a)


class ColourScale(unittest.TestCase):
    def test_default_follows_the_embedding_and_is_not_written(self):
        self.assertNotIn("--normalize_umap", argv(normalize="zscore"))

    def test_an_explicit_scale_is_passed_in_both_layouts(self):
        for layout in ("markers", "gates"):
            a = argv(layout=layout, normalize_umap="none")
            self.assertEqual(a[a.index("--normalize_umap") + 1], "none")


class GateSeparation(unittest.TestCase):
    def test_off_by_default(self):
        self.assertNotIn("--gate_separation", argv())

    def test_markers_layout_passes_it(self):
        a = argv(gate_separation=0.75)
        self.assertEqual(a[a.index("--gate_separation") + 1], "0.75")

    def test_slider_floats(self):
        self.assertNotIn("--gate_separation", argv(gate_separation=0.0))
        a = argv(gate_separation=1.0)
        self.assertEqual(float(a[a.index("--gate_separation") + 1]), 1.0)

    def test_gates_layout_does_not(self):
        self.assertNotIn("--gate_separation", argv(layout="gates", gate_separation=1))


class PointSize(unittest.TestCase):
    def test_auto_sends_nothing(self):
        a = argv()
        self.assertNotIn("--point_size", a)
        self.assertNotIn("--point_alpha", a)

    def test_a_chosen_size_is_sent(self):
        a = argv(point_size=0.1)
        self.assertEqual(a[a.index("--point_size") + 1], "0.1")


class Clip(unittest.TestCase):
    def test_off_by_default(self):
        self.assertNotIn("--clip_quantile", argv())

    def test_chosen_in_both_layouts(self):
        for layout in ("markers", "gates"):
            a = argv(layout=layout, clip_quantile="0.001")
            self.assertEqual(a[a.index("--clip_quantile") + 1], "0.001")


if __name__ == "__main__":
    unittest.main()
