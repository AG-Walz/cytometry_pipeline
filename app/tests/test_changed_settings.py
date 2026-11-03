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


if __name__ == "__main__":
    unittest.main()
