"""Unit tests for the v0.14 per-row tier classifier in bin/phase5_score.py."""
from __future__ import annotations

import sys
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "bin"))

from phase5_score import tier_for_row  # noqa: E402


class TestTierForRow(unittest.TestCase):
    def test_non_driver_is_untiered(self):
        self.assertEqual(tier_for_row("VUS", "approved"), "")
        self.assertEqual(tier_for_row("PASSENGER", "approved"), "")
        self.assertEqual(tier_for_row("OFF_TARGET", "phase2"), "")
        self.assertEqual(tier_for_row("", "approved"), "")

    def test_driver_approved_is_tier_1(self):
        self.assertEqual(tier_for_row("DRIVER", "approved"), "1")

    def test_driver_late_phase_is_tier_2(self):
        self.assertEqual(tier_for_row("DRIVER", "phase2"), "2")
        self.assertEqual(tier_for_row("DRIVER", "phase3"), "2")
        # case-insensitive
        self.assertEqual(tier_for_row("DRIVER", "Phase3"), "2")

    def test_driver_early_phase_is_tier_3(self):
        self.assertEqual(tier_for_row("DRIVER", "phase1"), "3")
        self.assertEqual(tier_for_row("DRIVER", "preclinical"), "3")

    def test_driver_unknown_phase_defaults_to_tier_3(self):
        # Conservative: unknown -> Tier 3 rather than untiered.
        self.assertEqual(tier_for_row("DRIVER", ""), "3")
        self.assertEqual(tier_for_row("DRIVER", "exploratory"), "3")


class TestPhase5OutputColumns(unittest.TestCase):
    """Black-box: a real phase 1-5 run on a toy DRIVER fixture must emit a
    'tier' column with the right value."""

    def test_toy_fgfr4_pipeline_writes_tier_column(self):
        import csv
        import subprocess
        import tempfile

        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            p1 = tmp_path / "p1.tsv"
            p2 = tmp_path / "p2.tsv"
            p3 = tmp_path / "p3.tsv"
            p4 = tmp_path / "p4.tsv"
            p5_tsv = tmp_path / "p5.tsv"
            p5_md = tmp_path / "report.md"

            for cmd in (
                ["python3", "bin/phase1_annotate.py",
                 "--vcf", "tests/data/toy_fgfr4_only.vcf",
                 "--cna", "assets/empty.cna.tsv",
                 "--fusion", "assets/empty.fusion.tsv",
                 "--targets-kb", "assets/targets_kb.tsv",
                 "--sample-id", "TOY_FGFR4",
                 "--out", str(p1)],
                ["python3", "bin/phase2_structure.py", "--in", str(p1), "--out", str(p2)],
                ["python3", "bin/phase3_dependency.py", "--in", str(p2),
                 "--depmap", "assets/depmap_rms_summary.tsv",
                 "--subtype", "FN", "--out", str(p3)],
                ["python3", "bin/phase4_drugs.py", "--in", str(p3),
                 "--drug-map", "assets/drug_target_map.tsv",
                 "--out", str(p4)],
                ["python3", "bin/phase5_score.py", "--in", str(p4),
                 "--vcf", "tests/data/toy_fgfr4_only.vcf",
                 "--sample-id", "TOY_FGFR4", "--subtype", "FN",
                 "--out-tsv", str(p5_tsv), "--out-md", str(p5_md)],
            ):
                subprocess.run(cmd, check=True, cwd=REPO_ROOT)

            with p5_tsv.open() as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                rows = list(reader)

            self.assertGreater(len(rows), 0)
            self.assertIn("tier", rows[0])

            driver_rows = [r for r in rows if r.get("call") == "DRIVER"]
            self.assertGreater(len(driver_rows), 0,
                               "FGFR4 hotspot fixture should produce DRIVER rows")
            for r in driver_rows:
                self.assertIn(r["tier"], ("1", "2", "3"),
                              f"DRIVER row got unexpected tier {r['tier']!r}")

            non_driver = [r for r in rows if r.get("call") != "DRIVER"]
            for r in non_driver:
                self.assertEqual(r["tier"], "",
                                 f"non-DRIVER row got tier {r['tier']!r}")

            md = p5_md.read_text()
            self.assertIn("Therapeutic tier summary", md)
            # Tier 1 expected because FGFR4 V550L -> erdafitinib (approved).
            self.assertIn("Tier 1", md)


if __name__ == "__main__":
    unittest.main()
