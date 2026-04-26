"""Unit tests for bin/cohort_visualize.py."""
from __future__ import annotations

import sys
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "bin"))

import cohort_visualize  # noqa: E402

FIXTURE_DIR = REPO_ROOT / "tests" / "data" / "cohort_viz"


class TestModuleSurface(unittest.TestCase):
    def test_constants_present(self):
        self.assertEqual(cohort_visualize.DRUGGABILITY_THRESHOLD, 0.10)
        self.assertEqual(len(cohort_visualize.HEATMAP_RAMP), 5)
        self.assertEqual(len(cohort_visualize.HEATMAP_BINS), 5)
        self.assertEqual(cohort_visualize.PER_SAMPLE_MAX_N, 100)

    def test_subtype_palette(self):
        self.assertEqual(set(cohort_visualize.SUBTYPE_PALETTE),
                         {"FN", "FP", "ALL"})


if __name__ == "__main__":
    unittest.main()
