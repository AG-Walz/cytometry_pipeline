"""The settings record in the results view. Run from the repo root:
python -m unittest app.tests.test_results_settings"""

import json
import tempfile
import unittest
from pathlib import Path

from app.services import results_scan


class SettingsRecord(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.run = Path(self.tmp.name)

    def test_record_is_read_and_offered_for_download(self):
        run_dir = self.run / "_run"
        run_dir.mkdir()
        (run_dir / "settings.json").write_text(json.dumps(
            {"status": "failed", "error": "population not found",
             "git": {"describe": "abc123 x", "dirty": True}}))
        (run_dir / "command.sh").write_text("Rscript s.R --n_neighbors='15'\n")
        r = results_scan.scan(self.run)
        self.assertEqual(r.settings["status"], "failed")
        self.assertIn("--n_neighbors='15'", r.command)
        self.assertEqual(sorted(p.name for p in r.settings_files),
                         ["command.sh", "settings.json"])
        # _run/ is not scanned for plots or data tables.
        self.assertEqual(r.data_files, [])

    def test_run_without_record_has_none(self):
        r = results_scan.scan(self.run)
        self.assertEqual((r.settings, r.command, r.settings_files), ({}, "", []))

    def test_section_is_built_for_a_failed_run_without_data(self):
        from app.callbacks import results
        (self.run / "_run").mkdir()
        (self.run / "_run" / "command.sh").write_text("Rscript s.R\n")
        (self.run / "_run" / "settings.json").write_text(
            json.dumps({"status": "failed", "error": "boom"}))
        self.assertTrue(results._data_section(results_scan.scan(self.run)))


class EmbeddingCheck(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.run = Path(self.tmp.name)

    def test_suspect_embedding_is_read(self):
        (self.run / "run_info.json").write_text(json.dumps(
            {"embedding": {"suspect": True, "reasons": ["The map is 550 units wide"]}}))
        r = results_scan.scan(self.run)
        self.assertTrue(r.embedding["suspect"])

    def test_no_run_info_means_no_verdict(self):
        self.assertEqual(results_scan.scan(self.run).embedding, {})


if __name__ == "__main__":
    unittest.main()
