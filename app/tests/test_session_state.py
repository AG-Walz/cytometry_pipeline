"""Form start and remembered sessions. Run from the repo root:
python -m unittest app.tests.test_session_state"""

import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from app import config
from app.services import session_state


class FormStart(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.dir.cleanup)
        self.file = Path(self.dir.name) / "session.json"
        for target, value in (("SESSION_FILE", self.file),
                              ("CACHE_DIR", Path(self.dir.name))):
            p = mock.patch.object(session_state, target, value)
            p.start()
            self.addCleanup(p.stop)

    def write(self, payload):
        self.file.write_text(json.dumps(payload), encoding="utf-8")

    def test_fresh_session_starts_on_marker_values(self):
        state = session_state.load()
        self.assertEqual(state["layout"], "markers")
        self.assertEqual(state["min_dist"], 0.1)

    def test_old_remembered_gate_default_moves_to_marker_values(self):
        self.write({"layout": "gates", "min_dist": 0.25, "point_size": 0.3})
        state = session_state.load()
        self.assertEqual(state["layout"], "markers")
        self.assertEqual(state["min_dist"], 0.1)
        self.assertEqual(state["point_size"], 0.3)   # other choices kept

    def test_gates_chosen_after_the_change_is_kept(self):
        session_state.save({"layout": "gates", "min_dist": 0.25})
        self.assertEqual(session_state.load()["layout"], "gates")

    def test_old_point_defaults_become_automatic(self):
        self.write({"layout": "markers", "form_start_version": 2,
                    "point_size": 0.4, "point_alpha": 0.6})
        state = session_state.load()
        self.assertIsNone(state["point_size"])
        self.assertIsNone(state["point_alpha"])
        self.assertEqual(state["layout"], "markers")   # v2 choice kept

    def test_a_chosen_point_size_is_kept(self):
        self.write({"layout": "markers", "form_start_version": 2, "point_size": 0.1})
        self.assertEqual(session_state.load()["point_size"], 0.1)

    def test_gates_is_not_labelled_recommended(self):
        labels = dict(config.LAYOUTS)
        self.assertNotIn("recommend", labels["gates"].lower())
        self.assertEqual(config.LAYOUTS[0][0], "markers")


if __name__ == "__main__":
    unittest.main()
