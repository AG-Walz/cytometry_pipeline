"""Runner regressions. Run from the repo root: python -m unittest app.tests.test_runner"""

import shutil
import tempfile
import threading
import time
import unittest

from app.services import detached, runner


class RefreshAfterNewRunStarted(unittest.TestCase):
    """A poll that began reading an old detached run's log must not be applied
    to a run that started while it was reading."""

    def setUp(self):
        self.out = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.out, ignore_errors=True)
        self.run = runner.Runner()
        ok, _ = self.run.start_detached(["/bin/sh", "-c", "echo hello"],
                                        cwd="/", output_dir=self.out,
                                        kind="run")
        self.assertTrue(ok)
        for _ in range(100):
            if self.run.snapshot()["status"] != runner.RUNNING:
                break
            time.sleep(0.05)
        self.assertEqual(self.run.snapshot()["status"], runner.DONE)

        real = detached.read_log
        def slow(*a, **k):
            time.sleep(0.3)
            return real(*a, **k)
        detached.read_log = slow
        self.addCleanup(setattr, detached, "read_log", real)

    def _start_inspect_during_poll(self):
        poll = threading.Thread(target=self.run.snapshot)
        poll.start()
        time.sleep(0.1)
        ok, _ = self.run.start(["/bin/sh", "-c", "sleep 1; echo 'Load report: x'"],
                               cwd="/", kind="inspect")
        self.assertTrue(ok)
        poll.join()
        return self.run.snapshot()

    def test_inspect_keeps_running_when_old_dir_exists(self):
        snap = self._start_inspect_during_poll()
        self.assertEqual(snap["status"], runner.RUNNING)
        self.assertNotIn("hello", snap["log"])

    def test_inspect_not_failed_when_old_dir_was_deleted(self):
        shutil.rmtree(self.out)
        snap = self._start_inspect_during_poll()
        self.assertEqual(snap["status"], runner.RUNNING)
        self.assertIsNone(snap["returncode"])
        time.sleep(1.5)
        self.assertEqual(self.run.snapshot()["status"], runner.DONE)


if __name__ == "__main__":
    unittest.main()
