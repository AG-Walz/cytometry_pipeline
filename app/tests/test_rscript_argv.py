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


if __name__ == "__main__":
    unittest.main()
