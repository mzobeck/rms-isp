"""Unit tests for bin/opentargets_client.py."""
from __future__ import annotations

import io
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "bin"))

import opentargets_client as ot  # noqa: E402


def _disease_payload(name: str, count: int, gene_scores: dict[str, float]) -> dict:
    return {
        "data": {
            "disease": {
                "id": "EFO_TEST",
                "name": name,
                "associatedTargets": {
                    "count": count,
                    "rows": [
                        {"score": s, "target": {"id": f"ENSG_{g}",
                                                 "approvedSymbol": g}}
                        for g, s in gene_scores.items()
                    ],
                },
            }
        }
    }


class TestFetchDiseaseTargetScores(unittest.TestCase):
    def test_parses_gene_scores(self):
        with tempfile.TemporaryDirectory() as tmp:
            payload = _disease_payload("rhabdomyosarcoma", 3, {
                "FGFR4": 0.6, "TP53": 0.7, "CDK4": 0.4,
            })
            with mock.patch("opentargets_client.urllib.request.urlopen") as m:
                m.return_value.__enter__.return_value = io.BytesIO(
                    json.dumps(payload).encode())
                summary = ot.fetch_disease_target_scores(
                    "EFO_0002918", Path(tmp))
            self.assertEqual(summary["disease_name"], "rhabdomyosarcoma")
            self.assertEqual(summary["count"], 3)
            self.assertAlmostEqual(summary["gene_scores"]["FGFR4"], 0.6)
            self.assertAlmostEqual(summary["gene_scores"]["TP53"], 0.7)

    def test_cache_hit_avoids_network(self):
        with tempfile.TemporaryDirectory() as tmp:
            cpath = Path(tmp) / "EFO_0002918.json"
            cpath.write_text(json.dumps({
                "disease_name": "rms",
                "efo_id": "EFO_0002918",
                "count": 1,
                "gene_scores": {"FGFR4": 0.5},
            }))
            with mock.patch("opentargets_client.urllib.request.urlopen") as m:
                summary = ot.fetch_disease_target_scores(
                    "EFO_0002918", Path(tmp))
            self.assertEqual(m.call_count, 0)
            self.assertEqual(summary["gene_scores"]["FGFR4"], 0.5)

    def test_null_disease_returns_none(self):
        with tempfile.TemporaryDirectory() as tmp:
            payload = {"data": {"disease": None}}
            with mock.patch("opentargets_client.urllib.request.urlopen") as m:
                m.return_value.__enter__.return_value = io.BytesIO(
                    json.dumps(payload).encode())
                summary = ot.fetch_disease_target_scores(
                    "EFO_BOGUS", Path(tmp))
            self.assertIsNone(summary)

    def test_network_failure_returns_none(self):
        import urllib.error
        with tempfile.TemporaryDirectory() as tmp, \
             mock.patch("opentargets_client.urllib.request.urlopen",
                        side_effect=urllib.error.URLError("down")):
            summary = ot.fetch_disease_target_scores("EFO_X", Path(tmp))
        self.assertIsNone(summary)


class TestLookupGeneDisease(unittest.TestCase):
    def test_known_gene_returns_score(self):
        with tempfile.TemporaryDirectory() as tmp:
            payload = _disease_payload("rms", 2, {"FGFR4": 0.6, "TP53": 0.7})
            with mock.patch("opentargets_client.urllib.request.urlopen") as m:
                m.return_value.__enter__.return_value = io.BytesIO(
                    json.dumps(payload).encode())
                result = ot.lookup_gene_disease(
                    "FGFR4", "EFO_0002918", Path(tmp))
            self.assertAlmostEqual(result["association_score"], 0.6)
            self.assertEqual(result["matched_disease_name"], "rms")

    def test_unknown_gene_returns_zero(self):
        with tempfile.TemporaryDirectory() as tmp:
            payload = _disease_payload("rms", 2, {"FGFR4": 0.6})
            with mock.patch("opentargets_client.urllib.request.urlopen") as m:
                m.return_value.__enter__.return_value = io.BytesIO(
                    json.dumps(payload).encode())
                result = ot.lookup_gene_disease(
                    "MADE_UP_GENE", "EFO_0002918", Path(tmp))
            self.assertEqual(result["association_score"], 0.0)
            # No matched disease name when score is 0 (gene not in OT's list).
            self.assertEqual(result["matched_disease_name"], "")

    def test_disease_lookup_failure_returns_none(self):
        with tempfile.TemporaryDirectory() as tmp:
            payload = {"data": {"disease": None}}
            with mock.patch("opentargets_client.urllib.request.urlopen") as m:
                m.return_value.__enter__.return_value = io.BytesIO(
                    json.dumps(payload).encode())
                result = ot.lookup_gene_disease(
                    "FGFR4", "EFO_BOGUS", Path(tmp))
            self.assertIsNone(result)

    def test_second_gene_uses_cache(self):
        """Two lookups for different genes against the same disease should
        result in exactly one network call (one disease query, cached)."""
        with tempfile.TemporaryDirectory() as tmp:
            payload = _disease_payload("rms", 2, {"FGFR4": 0.6, "TP53": 0.7})
            with mock.patch("opentargets_client.urllib.request.urlopen") as m:
                m.return_value.__enter__.return_value = io.BytesIO(
                    json.dumps(payload).encode())
                ot.lookup_gene_disease("FGFR4", "EFO_X", Path(tmp))
            with mock.patch("opentargets_client.urllib.request.urlopen") as m2:
                ot.lookup_gene_disease("TP53", "EFO_X", Path(tmp))
            self.assertEqual(m2.call_count, 0)


if __name__ == "__main__":
    unittest.main()
