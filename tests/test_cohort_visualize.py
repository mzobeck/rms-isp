"""Unit tests for bin/cohort_visualize.py."""
from __future__ import annotations

import csv
import shutil
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "bin"))

import cohort_visualize  # noqa: E402

FIXTURE_DIR = REPO_ROOT / "tests" / "data" / "cohort_viz"
SVG_NS = "{http://www.w3.org/2000/svg}"


class TestModuleSurface(unittest.TestCase):
    def test_constants_present(self):
        self.assertEqual(cohort_visualize.DRUGGABILITY_THRESHOLD, 0.10)
        self.assertEqual(len(cohort_visualize.HEATMAP_RAMP), 5)
        self.assertEqual(len(cohort_visualize.HEATMAP_BINS), 5)
        self.assertEqual(cohort_visualize.PER_SAMPLE_MAX_N, 100)

    def test_subtype_palette(self):
        self.assertEqual(set(cohort_visualize.SUBTYPE_PALETTE),
                         {"FN", "FP", "ALL"})


class TestAggregation(unittest.TestCase):
    def setUp(self):
        self.cohort_tsv = FIXTURE_DIR / "cohort_summary.tsv"
        self.target_rt_dir = FIXTURE_DIR

    def test_aggregate_returns_max_confidence_per_sample_gene(self):
        rows = cohort_visualize.aggregate_gene_matrix(
            self.cohort_tsv, self.target_rt_dir
        )
        a_fgfr4 = [r for r in rows if r["sample_id"] == "SAMP_A" and r["gene"] == "FGFR4"]
        self.assertEqual(len(a_fgfr4), 1)
        self.assertAlmostEqual(float(a_fgfr4[0]["max_confidence"]), 0.812)
        self.assertEqual(a_fgfr4[0]["subtype"], "FN")
        a_nras = [r for r in rows if r["sample_id"] == "SAMP_A" and r["gene"] == "NRAS"]
        self.assertEqual(a_nras, [])

    def test_aggregate_sort_is_deterministic(self):
        rows1 = cohort_visualize.aggregate_gene_matrix(
            self.cohort_tsv, self.target_rt_dir
        )
        rows2 = cohort_visualize.aggregate_gene_matrix(
            self.cohort_tsv, self.target_rt_dir
        )
        self.assertEqual(rows1, rows2)

    def test_write_gene_matrix_tsv_roundtrip(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "matrix.tsv"
            rows = cohort_visualize.aggregate_gene_matrix(
                self.cohort_tsv, self.target_rt_dir
            )
            cohort_visualize.write_gene_matrix_tsv(rows, out)
            with out.open() as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                read_rows = list(reader)
            self.assertEqual(len(read_rows), len(rows))
            self.assertEqual(set(read_rows[0]),
                             {"sample_id", "study", "subtype", "gene", "max_confidence"})


class TestMechanismChart(unittest.TestCase):
    def _bars(self, svg: str) -> list[ET.Element]:
        root = ET.fromstring(svg)
        return [r for r in root.iter(f"{SVG_NS}rect") if r.get("class") == "bar"]

    def test_one_bar_per_mechanism(self):
        svg = cohort_visualize.render_mechanism_chart(
            {"MEK_inhibitor": 5, "FGFR_inhibitor": 3, "CDK4/6_inhibitor": 2},
            cohort_size=10,
        )
        self.assertEqual(len(self._bars(svg)), 3)

    def test_long_tail_bins(self):
        svg = cohort_visualize.render_mechanism_chart(
            {"big": 50, "med": 10, "tiny1": 1, "tiny2": 1, "tiny3": 1},
            cohort_size=200,
        )
        bars = self._bars(svg)
        self.assertEqual(len(bars), 3)
        labels = [b.get("data-label") for b in bars]
        self.assertIn("other (3 items)", labels)

    def test_long_tail_omitted_when_empty(self):
        svg = cohort_visualize.render_mechanism_chart(
            {"big": 5, "med": 3}, cohort_size=10,
        )
        bars = self._bars(svg)
        self.assertEqual(len(bars), 2)
        self.assertNotIn("other", "".join(b.get("data-label", "") for b in bars))

    def test_output_is_valid_xml(self):
        svg = cohort_visualize.render_mechanism_chart({"x": 1}, cohort_size=1)
        root = ET.fromstring(svg)
        self.assertTrue(root.tag.endswith("svg"))


class TestDruggabilityChart(unittest.TestCase):
    def setUp(self):
        self.targets_kb = REPO_ROOT / "assets" / "targets_kb.tsv"

    def test_load_target_genes_returns_full_list(self):
        genes = cohort_visualize._load_target_genes(self.targets_kb)
        self.assertEqual(len(genes), 21)
        self.assertNotIn("gene", genes)
        self.assertFalse(any(g.startswith("#") for g in genes))

    def test_compute_druggability_matrix(self):
        agg_rows = [
            {"sample_id": "SAMP_A", "subtype": "FN", "gene": "FGFR4", "max_confidence": "0.812"},
            {"sample_id": "SAMP_B", "subtype": "FP", "gene": "CDK4",  "max_confidence": "0.701"},
            {"sample_id": "SAMP_B", "subtype": "FP", "gene": "FGFR4", "max_confidence": "0.42"},
            {"sample_id": "SAMP_C", "subtype": "FN", "gene": "NRAS",  "max_confidence": "0.654"},
        ]
        cohort_meta = [
            {"sample_id": "SAMP_A", "subtype": "FN"},
            {"sample_id": "SAMP_B", "subtype": "FP"},
            {"sample_id": "SAMP_C", "subtype": "FN"},
        ]
        m = cohort_visualize.compute_druggability_matrix(
            agg_rows, cohort_meta, genes=["FGFR4", "CDK4", "NRAS"],
        )
        self.assertAlmostEqual(m[("FGFR4", "FN")], 0.5)
        self.assertAlmostEqual(m[("FGFR4", "FP")], 1.0)
        self.assertAlmostEqual(m[("FGFR4", "ALL")], 0.0)
        self.assertAlmostEqual(m[("FGFR4", "whole_cohort")], 2 / 3)
        self.assertAlmostEqual(m[("CDK4", "FN")], 0.0)
        self.assertAlmostEqual(m[("CDK4", "FP")], 1.0)

    def test_color_bins(self):
        self.assertEqual(cohort_visualize._color_for_value(0.0),
                         cohort_visualize.EMPTY_CELL)
        self.assertEqual(cohort_visualize._color_for_value(0.05),
                         cohort_visualize.EMPTY_CELL)
        self.assertEqual(cohort_visualize._color_for_value(0.10), "#eff3ff")
        self.assertEqual(cohort_visualize._color_for_value(0.42), "#6baed6")
        self.assertEqual(cohort_visualize._color_for_value(0.99), "#08519c")
        self.assertEqual(cohort_visualize._color_for_value(1.0), "#08519c")

    def test_render_druggability_has_84_cells(self):
        agg_rows = [{"sample_id": "S1", "subtype": "FN", "gene": "FGFR4",
                     "max_confidence": "0.5"}]
        cohort_meta = [{"sample_id": "S1", "subtype": "FN"}]
        genes = cohort_visualize._load_target_genes(self.targets_kb)
        m = cohort_visualize.compute_druggability_matrix(
            agg_rows, cohort_meta, genes=genes,
        )
        svg = cohort_visualize.render_druggability_chart(
            m, genes, ["FN", "FP", "ALL", "whole_cohort"]
        )
        root = ET.fromstring(svg)
        cells = [r for r in root.iter(f"{SVG_NS}rect") if r.get("class") == "cell"]
        self.assertEqual(len(cells), 21 * 4)


class TestPerSampleHeatmap(unittest.TestCase):
    def _agg_and_meta(self, n_samples: int) -> tuple[list[dict], list[dict]]:
        agg = []
        meta = []
        for i in range(n_samples):
            sid = f"S{i:04d}"
            sub = "FN" if i % 2 == 0 else "FP"
            meta.append({"sample_id": sid, "subtype": sub})
            agg.append({"sample_id": sid, "subtype": sub,
                        "gene": "FGFR4", "max_confidence": "0.5"})
        return agg, meta

    def test_renders_for_small_cohort(self):
        agg, meta = self._agg_and_meta(3)
        svg = cohort_visualize.render_per_sample_heatmap(
            agg, meta, genes=["FGFR4", "CDK4"],
        )
        self.assertIsInstance(svg, str)
        root = ET.fromstring(svg)
        cells = [r for r in root.iter(f"{SVG_NS}rect") if r.get("class") == "cell"]
        self.assertEqual(len(cells), 2 * 3)

    def test_returns_none_above_threshold(self):
        agg, meta = self._agg_and_meta(150)
        result = cohort_visualize.render_per_sample_heatmap(
            agg, meta, genes=["FGFR4"],
        )
        self.assertIsNone(result)

    def test_subtype_stripe_present(self):
        agg, meta = self._agg_and_meta(4)
        svg = cohort_visualize.render_per_sample_heatmap(
            agg, meta, genes=["FGFR4"],
        )
        root = ET.fromstring(svg)
        stripes = [r for r in root.iter(f"{SVG_NS}rect")
                   if r.get("class") == "subtype-stripe"]
        self.assertEqual(len(stripes), 4)
        fills = [s.get("fill") for s in stripes]
        self.assertEqual(fills.count(cohort_visualize.SUBTYPE_PALETTE["FN"]), 2)
        self.assertEqual(fills.count(cohort_visualize.SUBTYPE_PALETTE["FP"]), 2)


class TestMain(unittest.TestCase):
    def test_main_emits_three_files_and_returns_status(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            shutil.copy(FIXTURE_DIR / "cohort_summary.tsv",
                        tmp_path / "cohort_summary.tsv")
            for sid in ("SAMP_A", "SAMP_B", "SAMP_C"):
                (tmp_path / sid).mkdir()
                shutil.copy(FIXTURE_DIR / sid / "p5.tsv",
                            tmp_path / sid / "p5.tsv")
            out_dir = tmp_path / "out"
            status = cohort_visualize.main(
                cohort_tsv=tmp_path / "cohort_summary.tsv",
                target_rt_dir=tmp_path,
                out_dir=out_dir,
                targets_kb=REPO_ROOT / "assets" / "targets_kb.tsv",
            )
            self.assertEqual(status["n_samples"], 3)
            self.assertTrue(status["mechanisms"].exists())
            self.assertTrue(status["druggability"].exists())
            self.assertIsNotNone(status["per_sample"])
            self.assertTrue(status["per_sample"].exists())
            self.assertTrue((out_dir / "cohort_gene_matrix.tsv").exists())

    def test_main_suppresses_per_sample_at_large_n(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            cohort_path = tmp_path / "cohort_summary.tsv"
            cols = ["sample_id", "study", "subtype", "pax_fusion", "histology",
                    "n_muts", "n_cnas", "n_fusions", "top_gene", "top_event_type",
                    "top_event", "top_drug", "top_mechanism", "top_confidence",
                    "top_call"]
            with cohort_path.open("w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
                w.writeheader()
                for i in range(150):
                    sid = f"S{i:04d}"
                    w.writerow({"sample_id": sid, "study": "x",
                                "subtype": "FN" if i % 2 else "FP",
                                "pax_fusion": "", "histology": "",
                                "n_muts": "1", "n_cnas": "0", "n_fusions": "0",
                                "top_gene": "FGFR4", "top_event_type": "snv",
                                "top_event": "V550L", "top_drug": "erdafitinib",
                                "top_mechanism": "pan_FGFR_inhibitor",
                                "top_confidence": "0.5",
                                "top_call": "DRIVER"})
                    p5_dir = tmp_path / sid
                    p5_dir.mkdir()
                    with (p5_dir / "p5.tsv").open("w", newline="") as p5fh:
                        p5w = csv.DictWriter(
                            p5fh,
                            fieldnames=["sample_id", "gene", "confidence"],
                            delimiter="\t",
                        )
                        p5w.writeheader()
                        p5w.writerow({"sample_id": sid, "gene": "FGFR4",
                                      "confidence": "0.5"})
            out_dir = tmp_path / "out"
            status = cohort_visualize.main(
                cohort_tsv=cohort_path,
                target_rt_dir=tmp_path,
                out_dir=out_dir,
                targets_kb=REPO_ROOT / "assets" / "targets_kb.tsv",
            )
            self.assertEqual(status["n_samples"], 150)
            self.assertIsNone(status["per_sample"])
            self.assertFalse((out_dir / "cohort_per_sample.svg").exists())


if __name__ == "__main__":
    unittest.main()
