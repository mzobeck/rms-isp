"""Unit tests for bin/cohort_portfolio.py."""
from __future__ import annotations

import csv
import sys
import tempfile
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "bin"))

import cohort_portfolio  # noqa: E402


def _write_p5(path: Path, rows: list[dict]) -> None:
    cols = ["sample_id", "gene", "call", "drug", "drug_mechanism",
            "drug_max_phase", "drug_pediatric_evidence", "confidence",
            "hgvsp_short", "consequence", "tier"]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in cols})


def _write_cohort_summary(path: Path, samples: list[dict]) -> None:
    cols = ["sample_id", "subtype"]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for s in samples:
            w.writerow({k: s.get(k, "") for k in cols})


class TestPortfolioTier(unittest.TestCase):
    def test_approved_above_5pct_is_tier_1(self):
        self.assertEqual(cohort_portfolio.portfolio_tier("approved", 0.10), "1")
        self.assertEqual(cohort_portfolio.portfolio_tier("approved", 0.05), "1")

    def test_approved_3_to_5pct_is_tier_2(self):
        self.assertEqual(cohort_portfolio.portfolio_tier("approved", 0.04), "2")
        self.assertEqual(cohort_portfolio.portfolio_tier("approved", 0.03), "2")

    def test_approved_below_3pct_is_tier_3(self):
        self.assertEqual(cohort_portfolio.portfolio_tier("approved", 0.02), "3")
        self.assertEqual(cohort_portfolio.portfolio_tier("approved", 0.0), "3")

    def test_late_phase_above_3pct_is_tier_2(self):
        self.assertEqual(cohort_portfolio.portfolio_tier("phase2", 0.10), "2")
        self.assertEqual(cohort_portfolio.portfolio_tier("phase3", 0.05), "2")

    def test_late_phase_below_3pct_is_tier_3(self):
        self.assertEqual(cohort_portfolio.portfolio_tier("phase2", 0.02), "3")

    def test_early_phase_always_tier_3(self):
        self.assertEqual(cohort_portfolio.portfolio_tier("phase1", 0.20), "3")
        self.assertEqual(cohort_portfolio.portfolio_tier("preclinical", 0.50), "3")

    def test_unknown_phase_is_unranked(self):
        self.assertEqual(cohort_portfolio.portfolio_tier("", 0.10), "")
        self.assertEqual(cohort_portfolio.portfolio_tier("exploratory", 0.10), "")


class TestPrevalenceMath(unittest.TestCase):
    def test_subtype_sample_sets(self):
        meta = [
            {"sample_id": "S1", "subtype": "FN"},
            {"sample_id": "S2", "subtype": "FN"},
            {"sample_id": "S3", "subtype": "FP"},
            {"sample_id": "S4", "subtype": "ALL"},
            {"sample_id": "S5", "subtype": "unknown"},
        ]
        sets = cohort_portfolio.subtype_sample_sets(meta)
        self.assertEqual(sets["FN"], {"S1", "S2"})
        self.assertEqual(sets["FP"], {"S3"})
        self.assertEqual(sets["ALL"], {"S4"})

    def test_gene_prevalence_division(self):
        drivers = {
            "FGFR4": {"FN": set(), "FP": {"S3"}, "ALL": set()},
            "TP53":  {"FN": {"S1", "S2"}, "FP": set(), "ALL": set()},
        }
        sizes = {"FN": 4, "FP": 2, "ALL": 1}
        prev = cohort_portfolio.gene_prevalence(drivers, sizes)
        self.assertAlmostEqual(prev["FGFR4"]["FP"], 0.5)
        self.assertAlmostEqual(prev["FGFR4"]["FN"], 0.0)
        self.assertAlmostEqual(prev["TP53"]["FN"], 0.5)


class TestEndToEndPortfolio(unittest.TestCase):
    """Build a synthetic 5-sample cohort and run the full pipeline."""

    def test_synthetic_portfolio(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            target_rt = tmp_path / "target_rt"
            target_rt.mkdir()

            samples = [
                {"sample_id": "FN1", "subtype": "FN"},
                {"sample_id": "FN2", "subtype": "FN"},
                {"sample_id": "FN3", "subtype": "FN"},
                {"sample_id": "FP1", "subtype": "FP"},
                {"sample_id": "FP2", "subtype": "FP"},
            ]
            cohort_tsv = tmp_path / "cohort_summary.tsv"
            _write_cohort_summary(cohort_tsv, samples)

            # FGFR4 hits in 2 of 3 FN samples -> 67% FN prevalence.
            # erdafitinib is approved -> Portfolio Tier 1.
            for sid in ("FN1", "FN2"):
                _write_p5(target_rt / sid / "p5.tsv", [
                    {"sample_id": sid, "gene": "FGFR4", "call": "DRIVER",
                     "drug": "erdafitinib", "drug_mechanism": "pan_FGFR_inhibitor",
                     "drug_max_phase": "approved",
                     "drug_pediatric_evidence": "yes_trial",
                     "confidence": "0.81", "hgvsp_short": "V550L"},
                ])

            # TP53 hit in 1 of 3 FN samples -> 33% FN prevalence.
            # WEE1i is phase 2 -> Portfolio Tier 2 (>= 3% prevalence).
            _write_p5(target_rt / "FN3" / "p5.tsv", [
                {"sample_id": "FN3", "gene": "TP53", "call": "DRIVER",
                 "drug": "adavosertib", "drug_mechanism": "WEE1_inhibitor",
                 "drug_max_phase": "phase2",
                 "drug_pediatric_evidence": "yes_trial",
                 "confidence": "0.55", "hgvsp_short": "R175H"},
            ])

            # CDK4 hit in 1 of 2 FP samples -> 50% FP prevalence.
            # palbociclib is approved -> Portfolio Tier 1.
            _write_p5(target_rt / "FP1" / "p5.tsv", [
                {"sample_id": "FP1", "gene": "CDK4", "call": "DRIVER",
                 "drug": "palbociclib", "drug_mechanism": "CDK4/6_inhibitor",
                 "drug_max_phase": "approved",
                 "drug_pediatric_evidence": "yes_trial",
                 "confidence": "0.72", "consequence": "amplification"},
            ])

            # FP2 has only a passenger - should not contribute to any tier.
            _write_p5(target_rt / "FP2" / "p5.tsv", [
                {"sample_id": "FP2", "gene": "TP53", "call": "PASSENGER",
                 "drug": "adavosertib", "drug_mechanism": "WEE1_inhibitor",
                 "drug_max_phase": "phase2",
                 "drug_pediatric_evidence": "yes_trial",
                 "confidence": "0.05"},
            ])

            out_path = tmp_path / "STS.md"
            status = cohort_portfolio.main(
                cohort_tsv=cohort_tsv,
                target_rt_dir=target_rt,
                out_path=out_path,
                pipeline_version="v0.15.0-test",
            )

            self.assertEqual(status["n_samples"], 5)
            self.assertEqual(status["n_portfolio_genes"], 3)  # FGFR4, TP53, CDK4
            self.assertEqual(status["tier_counts"]["1"], 2)   # FGFR4, CDK4
            self.assertEqual(status["tier_counts"]["2"], 1)   # TP53

            md = out_path.read_text()
            self.assertIn("FGFR4", md)
            self.assertIn("CDK4", md)
            self.assertIn("TP53", md)
            self.assertIn("Tier 1", md)
            self.assertIn("Tier 2", md)
            # Anchor table maps each sample to its highest-tier hit.
            self.assertIn("**FGFR4**", md)
            self.assertIn("`erdafitinib`", md)


class TestBestDrugSelection(unittest.TestCase):
    def test_higher_phase_wins(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            target_rt = tmp_path / "target_rt"
            target_rt.mkdir()
            samples = [{"sample_id": "S1", "subtype": "FN"}]
            cohort_tsv = tmp_path / "cohort_summary.tsv"
            _write_cohort_summary(cohort_tsv, samples)
            _write_p5(target_rt / "S1" / "p5.tsv", [
                {"sample_id": "S1", "gene": "X", "call": "DRIVER",
                 "drug": "preclinical_drug", "drug_max_phase": "preclinical",
                 "confidence": "0.5"},
                {"sample_id": "S1", "gene": "X", "call": "DRIVER",
                 "drug": "approved_drug", "drug_max_phase": "approved",
                 "confidence": "0.4"},
            ])
            drugs = cohort_portfolio.best_drug_per_gene(samples, target_rt)
            self.assertEqual(drugs["X"]["drug"], "approved_drug")
            self.assertEqual(drugs["X"]["max_phase"], "approved")

    def test_confidence_breaks_phase_tie(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            target_rt = tmp_path / "target_rt"
            target_rt.mkdir()
            samples = [{"sample_id": "S1", "subtype": "FN"}]
            cohort_tsv = tmp_path / "cohort_summary.tsv"
            _write_cohort_summary(cohort_tsv, samples)
            _write_p5(target_rt / "S1" / "p5.tsv", [
                {"sample_id": "S1", "gene": "X", "call": "DRIVER",
                 "drug": "low_conf", "drug_max_phase": "approved",
                 "confidence": "0.3"},
                {"sample_id": "S1", "gene": "X", "call": "DRIVER",
                 "drug": "high_conf", "drug_max_phase": "approved",
                 "confidence": "0.8"},
            ])
            drugs = cohort_portfolio.best_drug_per_gene(samples, target_rt)
            self.assertEqual(drugs["X"]["drug"], "high_conf")


if __name__ == "__main__":
    unittest.main()
