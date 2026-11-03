"""The version badge. Run from the repo root:
python -m unittest app.tests.test_version"""

import unittest
from unittest import mock

from app.services import version


def info(running, disk, remote, error=None):
    with mock.patch.object(version, "RUNNING", {"commit": running, "short": running[:7]}), \
         mock.patch.object(version, "_describe", return_value={"commit": disk}), \
         mock.patch.object(version, "_remote_main",
                           return_value={"commit": remote, "error": error}):
        return version.status()


class State(unittest.TestCase):
    def test_current(self):
        self.assertEqual(info("a" * 40, "a" * 40, "a" * 40)["state"], "current")

    def test_newer_checkout_on_disk_means_restart(self):
        self.assertEqual(info("a" * 40, "b" * 40, "b" * 40)["state"], "restart")

    def test_newer_github_means_behind(self):
        self.assertEqual(info("a" * 40, "a" * 40, "c" * 40)["state"], "behind")

    def test_github_unreachable_is_unknown_not_current(self):
        self.assertEqual(info("a" * 40, "a" * 40, None, "unreachable")["state"],
                         "unknown")

    def test_running_reads_this_checkout(self):
        self.assertEqual(len(version.RUNNING.get("commit", "")), 40)

    def test_badge_renders_for_every_state(self):
        from app.callbacks.version import badge
        for state in ("current", "restart", "behind", "unknown"):
            badge({"state": state, "running": {"short": "abc1234", "date": "2026-11-03 14:22"},
                   "on_disk": {"short": "def5678"}, "remote_commit": "f" * 40,
                   "remote_error": None})


if __name__ == "__main__":
    unittest.main()
