"""Unit tests for the noise-source filters in bin/fetch_dgidb.py."""
from __future__ import annotations

import sys
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "bin"))

import fetch_dgidb  # noqa: E402


class TestIsLowTrustSole(unittest.TestCase):
    def test_pure_low_trust_returns_true(self):
        self.assertTrue(fetch_dgidb.is_low_trust_sole(["ClearityFoundationClinicalTrial"]))

    def test_low_trust_plus_authoritative_returns_false(self):
        self.assertFalse(fetch_dgidb.is_low_trust_sole(
            ["ClearityFoundationClinicalTrial", "ChEMBL"]))

    def test_authoritative_only_returns_false(self):
        self.assertFalse(fetch_dgidb.is_low_trust_sole(["ChEMBL", "CIViC"]))

    def test_empty_returns_false(self):
        # Empty source list should not be classified as low-trust; the upstream
        # render_rows already handles missing sources separately.
        self.assertFalse(fetch_dgidb.is_low_trust_sole([]))


class TestRenderRowsFiltersGranisetron(unittest.TestCase):
    """The TP53 -> granisetron false-positive case should be dropped by the
    fetcher even though it has a typed mechanism."""

    def test_granisetron_dropped_by_low_trust_filter(self):
        payload = {
            "data": {
                "genes": {
                    "nodes": [
                        {
                            "name": "TP53",
                            "interactions": [
                                {
                                    "drug": {"name": "granisetron",
                                             "conceptId": "rxcui:26237",
                                             "approved": True},
                                    "interactionTypes": [{"type": "activator"}],
                                    "sources": [
                                        {"sourceDbName": "ClearityFoundationClinicalTrial"}
                                    ],
                                    "publications": [],
                                },
                                {
                                    "drug": {"name": "idasanutlin",
                                             "conceptId": "ncit:C99131",
                                             "approved": False},
                                    "interactionTypes": [{"type": "inhibitor"}],
                                    "sources": [
                                        {"sourceDbName": "ChEMBL"},
                                        {"sourceDbName": "CKB-CORE"},
                                    ],
                                    "publications": [{"pmid": "31483066"}],
                                },
                            ],
                        }
                    ]
                }
            }
        }
        rows, n_dropped_mech, n_dropped_low_trust = fetch_dgidb.render_rows(payload)
        drug_names = [r["drug"] for r in rows]
        self.assertNotIn("granisetron", drug_names)
        self.assertIn("idasanutlin", drug_names)
        self.assertEqual(n_dropped_mech, 0)
        self.assertEqual(n_dropped_low_trust, 1)

    def test_low_trust_plus_authoritative_kept(self):
        # If the low-trust source co-occurs with a real one, keep the row.
        payload = {
            "data": {
                "genes": {
                    "nodes": [
                        {
                            "name": "TP53",
                            "interactions": [
                                {
                                    "drug": {"name": "co_listed_drug",
                                             "conceptId": "x", "approved": False},
                                    "interactionTypes": [{"type": "inhibitor"}],
                                    "sources": [
                                        {"sourceDbName": "ClearityFoundationClinicalTrial"},
                                        {"sourceDbName": "ChEMBL"},
                                    ],
                                    "publications": [],
                                }
                            ],
                        }
                    ]
                }
            }
        }
        rows, _, n_low_trust = fetch_dgidb.render_rows(payload)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["drug"], "co_listed_drug")
        self.assertEqual(n_low_trust, 0)


if __name__ == "__main__":
    unittest.main()
