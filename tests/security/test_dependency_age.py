import importlib.util
from pathlib import Path
import unittest
from unittest.mock import patch
from datetime import datetime, timezone

spec = importlib.util.spec_from_file_location("age", Path(__file__).parents[2] / "scripts/check-fresh-deps.py")
age = importlib.util.module_from_spec(spec)
spec.loader.exec_module(age)


class DependencyAgeTests(unittest.TestCase):
    def test_missing_release_metadata_fails(self):
        with self.assertRaises(ValueError):
            age.release_time({"urls": []})

    def test_non_registry_source_rejected(self):
        with self.assertRaises(ValueError):
            age.registry_packages('[[package]]\nname="dep"\nversion="1"\nsource={git="https://example.com/repo"}')

    def test_network_failure_blocks(self):
        with patch.object(age.urllib.request, "urlopen", side_effect=TimeoutError):
            self.assertIn("cannot verify", age.check_package(("dep", "1"), datetime.now(timezone.utc)))

    def test_age_boundary(self):
        now = datetime(2026, 10, 7, tzinfo=timezone.utc)
        for timestamp, blocked in [("2026-09-30T00:00:00Z", False), ("2026-09-30T00:00:01Z", True)]:
            with patch.object(age.urllib.request, "urlopen"), patch.object(age.json, "load", return_value={"urls": [{"upload_time_iso_8601": timestamp}]}):
                self.assertEqual(age.check_package(("dep", "1"), now) is not None, blocked)
