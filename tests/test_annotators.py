"""Unit tests for bin/annotators/."""
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

from annotators import Variant, get_annotator  # noqa: E402
from annotators.base import VariantAnnotation, _Chain  # noqa: E402
from annotators.curated_vcf import CuratedVCFAnnotator, extract_hotspot  # noqa: E402
from annotators.oncokb import OncoKBAnnotator  # noqa: E402
from annotators.vep_rest import (  # noqa: E402
    BATCH_SIZE,
    VEPRestAnnotator,
    hgvsp_short,
    parse_vep_response,
    pick_canonical,
)


class TestFactory(unittest.TestCase):
    def test_get_curated(self):
        ann = get_annotator("curated")
        self.assertIsInstance(ann, CuratedVCFAnnotator)

    def test_get_vep_rest_with_cache_dir(self):
        with tempfile.TemporaryDirectory() as tmp:
            ann = get_annotator("vep_rest", cache_dir=Path(tmp))
            self.assertIsInstance(ann, VEPRestAnnotator)
            self.assertEqual(ann.cache_dir, Path(tmp))

    def test_unknown_raises(self):
        with self.assertRaises(ValueError):
            get_annotator("not_a_real_backend")

    def test_get_oncokb_uses_env_token(self):
        with tempfile.TemporaryDirectory() as tmp, \
             mock.patch.dict("os.environ", {"ONCOKB_TOKEN": "env-token"}, clear=False):
            ann = get_annotator("oncokb", oncokb_cache_dir=Path(tmp))
            self.assertIsInstance(ann, OncoKBAnnotator)
            self.assertEqual(ann.token, "env-token")

    def test_chain_returns_chain_instance(self):
        with tempfile.TemporaryDirectory() as tmp:
            ann = get_annotator(
                "vep_rest,oncokb",
                cache_dir=Path(tmp) / "vep",
                oncokb_cache_dir=Path(tmp) / "oncokb",
                oncokb_token="x",
            )
            self.assertIsInstance(ann, _Chain)
            self.assertEqual(len(ann.members), 2)
            self.assertIsInstance(ann.members[0], VEPRestAnnotator)
            self.assertIsInstance(ann.members[1], OncoKBAnnotator)

    def test_chain_passes_oncokb_tumor_type(self):
        with tempfile.TemporaryDirectory() as tmp:
            ann = get_annotator(
                "oncokb",
                oncokb_cache_dir=Path(tmp),
                oncokb_token="x",
                oncokb_tumor_type="ERMS",
            )
            self.assertEqual(ann.tumor_type, "ERMS")


class TestCuratedVCFAnnotator(unittest.TestCase):
    def test_pulls_gene_consequence_hgvsp_from_info(self):
        ann = CuratedVCFAnnotator()
        v = Variant(chrom="18", pos=70522538, ref="C", alt="T", info={
            "GENE": "FGFR4",
            "CONSEQUENCE": "missense",
            "NOTE": "kinase domain hotspot V550L per Shern 2014",
        })
        result = ann.annotate_batch([v])[0]
        self.assertEqual(result.gene, "FGFR4")
        self.assertEqual(result.consequence, "missense")
        self.assertEqual(result.hgvsp_short, "V550L")
        self.assertEqual(result.source, "curated")

    def test_empty_info_yields_empty_annotation(self):
        ann = CuratedVCFAnnotator()
        v = Variant(chrom="1", pos=1, ref="A", alt="T", info={})
        result = ann.annotate_batch([v])[0]
        self.assertEqual(result.gene, "")
        self.assertEqual(result.consequence, "")
        self.assertEqual(result.hgvsp_short, "")

    def test_extract_hotspot_isolates_pattern(self):
        self.assertEqual(extract_hotspot("R175H"), "R175H")
        self.assertEqual(extract_hotspot("hotspot R175H per literature"), "R175H")
        self.assertEqual(extract_hotspot("missense"), "")
        self.assertEqual(extract_hotspot(""), "")

    def test_batch_preserves_order(self):
        ann = CuratedVCFAnnotator()
        variants = [
            Variant("1", i, "A", "T", info={"GENE": f"GENE{i}"})
            for i in range(5)
        ]
        results = ann.annotate_batch(variants)
        for v, r in zip(variants, results):
            self.assertEqual(r.gene, v.info["GENE"])


class TestVEPRestHelpers(unittest.TestCase):
    def test_hgvsp_short_substitution(self):
        self.assertEqual(hgvsp_short("ENST00000379368.4:p.Val550Leu"), "V550L")
        self.assertEqual(hgvsp_short("p.Arg175His"), "R175H")
        self.assertEqual(hgvsp_short("p.Gln61Lys"), "Q61K")

    def test_hgvsp_short_returns_empty_for_unparseable(self):
        # Insertions, frameshifts, and other complex hgvsp formats.
        self.assertEqual(hgvsp_short(""), "")
        self.assertEqual(hgvsp_short("c.123A>G"), "")
        self.assertEqual(hgvsp_short("p.fakeval"), "")
        self.assertEqual(hgvsp_short("p.Val550Trpfs*23"), "")  # frameshift

    def test_pick_canonical_prefers_canonical_flag(self):
        tcs = [
            {"transcript_id": "ENST_a", "canonical": 0},
            {"transcript_id": "ENST_b", "canonical": 1},
            {"transcript_id": "ENST_c"},
        ]
        self.assertEqual(pick_canonical(tcs)["transcript_id"], "ENST_b")

    def test_pick_canonical_falls_back_to_first(self):
        tcs = [
            {"transcript_id": "ENST_a"},
            {"transcript_id": "ENST_b"},
        ]
        self.assertEqual(pick_canonical(tcs)["transcript_id"], "ENST_a")

    def test_pick_canonical_handles_empty(self):
        self.assertIsNone(pick_canonical([]))

    def test_parse_vep_response_extracts_canonical_fields(self):
        obj = {
            "most_severe_consequence": "missense_variant",
            "transcript_consequences": [
                {
                    "canonical": 1,
                    "gene_symbol": "FGFR4",
                    "consequence_terms": ["missense_variant"],
                    "hgvsp": "ENST00000379368.4:p.Val550Leu",
                },
            ],
        }
        gene, cons, hgvsp = parse_vep_response(obj)
        self.assertEqual(gene, "FGFR4")
        self.assertEqual(cons, "missense_variant")
        self.assertEqual(hgvsp, "V550L")

    def test_parse_vep_response_handles_no_transcripts(self):
        obj = {"most_severe_consequence": "intergenic_variant",
               "transcript_consequences": []}
        gene, cons, hgvsp = parse_vep_response(obj)
        self.assertEqual(gene, "")
        self.assertEqual(cons, "intergenic_variant")
        self.assertEqual(hgvsp, "")


def _canned_response(gene: str, cons: str, hgvsp: str) -> dict:
    return {
        "most_severe_consequence": cons,
        "transcript_consequences": [{
            "canonical": 1,
            "gene_symbol": gene,
            "consequence_terms": [cons],
            "hgvsp": hgvsp,
        }],
    }


class TestVEPRestAnnotator(unittest.TestCase):
    def test_batch_call_then_cache_hit(self):
        with tempfile.TemporaryDirectory() as tmp:
            ann = VEPRestAnnotator(cache_dir=Path(tmp))
            variants = [
                Variant("18", 70522538, "C", "T", info={}),
                Variant("17", 7675088, "G", "A", info={}),
            ]
            canned = [
                _canned_response("FGFR4", "missense_variant",
                                 "ENST00000379368.4:p.Val550Leu"),
                _canned_response("TP53", "missense_variant",
                                 "ENST00000269305.4:p.Arg175His"),
            ]
            response_body = json.dumps(canned).encode()

            with mock.patch("annotators.vep_rest.urllib.request.urlopen") as m_open:
                m_open.return_value.__enter__.return_value = io.BytesIO(response_body)
                results = ann.annotate_batch(variants)

            self.assertEqual(m_open.call_count, 1)
            self.assertEqual(len(results), 2)
            self.assertEqual(results[0].gene, "FGFR4")
            self.assertEqual(results[0].consequence, "missense_variant")
            self.assertEqual(results[0].hgvsp_short, "V550L")
            self.assertEqual(results[0].source, "vep_rest")
            self.assertEqual(results[1].gene, "TP53")
            self.assertEqual(results[1].hgvsp_short, "R175H")

            # Cache files written; second call should not hit the network.
            self.assertTrue((Path(tmp) / "18_70522538_C_T.json").exists())
            self.assertTrue((Path(tmp) / "17_7675088_G_A.json").exists())

            with mock.patch("annotators.vep_rest.urllib.request.urlopen") as m_open2:
                results2 = ann.annotate_batch(variants)
            self.assertEqual(m_open2.call_count, 0)
            self.assertEqual([r.gene for r in results], [r.gene for r in results2])

    def test_network_failure_degrades_gracefully(self):
        import urllib.error
        with tempfile.TemporaryDirectory() as tmp:
            ann = VEPRestAnnotator(cache_dir=Path(tmp))
            variants = [Variant("1", 100, "A", "T", info={})]
            with mock.patch("annotators.vep_rest.urllib.request.urlopen",
                            side_effect=urllib.error.URLError("network down")):
                results = ann.annotate_batch(variants)
            self.assertEqual(len(results), 1)
            self.assertEqual(results[0].consequence, "")
            self.assertEqual(results[0].source, "vep_rest")

    def test_non_dict_payload_degrades(self):
        with tempfile.TemporaryDirectory() as tmp:
            ann = VEPRestAnnotator(cache_dir=Path(tmp))
            variants = [Variant("1", 100, "A", "T", info={})]
            with mock.patch("annotators.vep_rest.urllib.request.urlopen") as m_open:
                m_open.return_value.__enter__.return_value = io.BytesIO(b'"not a list"')
                results = ann.annotate_batch(variants)
            self.assertEqual(results[0].consequence, "")

    def test_partial_cache_hit_only_fetches_misses(self):
        with tempfile.TemporaryDirectory() as tmp:
            ann = VEPRestAnnotator(cache_dir=Path(tmp))
            # Pre-populate one cache entry.
            cached = _canned_response("FGFR4", "missense_variant",
                                      "ENST00000379368.4:p.Val550Leu")
            (Path(tmp) / "18_70522538_C_T.json").write_text(json.dumps(cached))

            variants = [
                Variant("18", 70522538, "C", "T", info={}),  # cache hit
                Variant("1", 114713908, "G", "T", info={}),  # cache miss
            ]
            miss_response = json.dumps([
                _canned_response("NRAS", "missense_variant",
                                 "ENST00000369535.4:p.Gln61Lys"),
            ]).encode()
            with mock.patch("annotators.vep_rest.urllib.request.urlopen") as m_open:
                m_open.return_value.__enter__.return_value = io.BytesIO(miss_response)
                results = ann.annotate_batch(variants)
            self.assertEqual(m_open.call_count, 1)
            self.assertEqual(results[0].gene, "FGFR4")
            self.assertEqual(results[1].gene, "NRAS")

    def test_batch_size_constant_is_reasonable(self):
        # Sanity: Ensembl documents 200 as the batch limit.
        self.assertEqual(BATCH_SIZE, 200)


def _oncokb_response(oncogenic: str, effect: str) -> dict:
    return {"oncogenic": oncogenic, "mutationEffect": {"knownEffect": effect}}


class TestOncoKBAnnotator(unittest.TestCase):
    def _prior(self) -> tuple[list[Variant], list[VariantAnnotation]]:
        variants = [
            Variant("18", 70522538, "C", "T"),
            Variant("17", 7675088, "G", "A"),
            Variant("3", 178917478, "T", "C"),
        ]
        prior = [
            VariantAnnotation(gene="FGFR4", consequence="missense_variant",
                              hgvsp_short="V550L", source="vep_rest"),
            VariantAnnotation(gene="TP53", consequence="missense_variant",
                              hgvsp_short="R175H", source="vep_rest"),
            # Empty hgvsp -> ineligible, should pass through untouched.
            VariantAnnotation(gene="PIK3CA", consequence="intron_variant",
                              hgvsp_short="", source="vep_rest"),
        ]
        return variants, prior

    def test_no_token_passes_through(self):
        with tempfile.TemporaryDirectory() as tmp:
            ann = OncoKBAnnotator(cache_dir=Path(tmp), token=None)
            ann.token = ""  # explicit override of any env-leaked value
            variants, prior = self._prior()
            with mock.patch.dict("os.environ", {}, clear=True):
                result = ann.annotate_with_prior(variants, prior)
            self.assertEqual(result, prior)

    def test_with_token_lifts_oncogenic(self):
        with tempfile.TemporaryDirectory() as tmp:
            ann = OncoKBAnnotator(cache_dir=Path(tmp), token="abc",
                                   tumor_type="RMS")
            variants, prior = self._prior()
            canned = [
                _oncokb_response("Likely Oncogenic", "Gain-of-function"),
                _oncokb_response("Oncogenic", "Loss-of-function"),
            ]
            response_body = json.dumps(canned).encode()
            with mock.patch("annotators.oncokb.urllib.request.urlopen") as m_open:
                m_open.return_value.__enter__.return_value = io.BytesIO(response_body)
                result = ann.annotate_with_prior(variants, prior)

            self.assertEqual(m_open.call_count, 1)
            # Eligible entries got augmented.
            self.assertEqual(result[0].oncogenic, "Likely Oncogenic")
            self.assertEqual(result[0].mutation_effect, "Gain-of-function")
            self.assertEqual(result[0].source, "vep_rest+oncokb")
            self.assertEqual(result[1].oncogenic, "Oncogenic")
            # Ineligible (empty hgvsp_short) was NOT sent to OncoKB and stays put.
            self.assertEqual(result[2], prior[2])

    def test_cache_hit_avoids_network(self):
        with tempfile.TemporaryDirectory() as tmp:
            ann = OncoKBAnnotator(cache_dir=Path(tmp), token="abc",
                                   tumor_type="RMS")
            cached_path = Path(tmp) / "FGFR4__V550L__RMS.json"
            cached_path.write_text(json.dumps(
                _oncokb_response("Likely Oncogenic", "Gain-of-function")))

            variants = [Variant("18", 70522538, "C", "T")]
            prior = [VariantAnnotation(gene="FGFR4", consequence="missense_variant",
                                        hgvsp_short="V550L", source="vep_rest")]
            with mock.patch("annotators.oncokb.urllib.request.urlopen") as m_open:
                result = ann.annotate_with_prior(variants, prior)

            self.assertEqual(m_open.call_count, 0)
            self.assertEqual(result[0].oncogenic, "Likely Oncogenic")

    def test_network_failure_passes_prior_through(self):
        import urllib.error
        with tempfile.TemporaryDirectory() as tmp:
            ann = OncoKBAnnotator(cache_dir=Path(tmp), token="abc")
            variants = [Variant("18", 70522538, "C", "T")]
            prior = [VariantAnnotation(gene="FGFR4", consequence="missense_variant",
                                        hgvsp_short="V550L", source="vep_rest")]
            with mock.patch("annotators.oncokb.urllib.request.urlopen",
                            side_effect=urllib.error.URLError("oncokb down")):
                result = ann.annotate_with_prior(variants, prior)
            self.assertEqual(result, prior)


class TestChainAnnotators(unittest.TestCase):
    def test_chain_runs_first_pass_then_augments(self):
        from annotators.base import _Chain

        class FakeFirst:
            def annotate_batch(self, variants):
                return [VariantAnnotation(gene="FGFR4", consequence="missense_variant",
                                           hgvsp_short="V550L", source="fake_first")
                        for _ in variants]

        class FakeSecond:
            def annotate_with_prior(self, variants, prior):
                return [VariantAnnotation(
                    gene=p.gene, consequence=p.consequence, hgvsp_short=p.hgvsp_short,
                    source=p.source + "+fake_second",
                    oncogenic="Oncogenic", mutation_effect="Gain-of-function",
                ) for p in prior]

        chain = _Chain([FakeFirst(), FakeSecond()])
        result = chain.annotate_batch([Variant("1", 1, "A", "T")])
        self.assertEqual(result[0].oncogenic, "Oncogenic")
        self.assertEqual(result[0].source, "fake_first+fake_second")

    def test_classify_snv_lift_on_oncokb_oncogenic(self):
        # If OncoKB calls a non-hotspot missense Oncogenic, classify_snv lifts
        # to DRIVER. Critical to v0.13 behavior.
        sys.path.insert(0, str(REPO_ROOT / "bin"))
        from phase1_annotate import classify_snv

        # FGFR4 has hotspots in the targets KB.
        gene_kb = {"hotspots": {"V550L"}, "lof_target": False,
                   "amp_target": False, "del_target": False, "role": "oncogene"}

        # Variant not on hotspot list, but OncoKB-Oncogenic -> DRIVER.
        call, reason, score = classify_snv(
            "missense_variant", gene_kb, "Y367C", oncogenic="Likely Oncogenic",
        )
        self.assertEqual(call, "DRIVER")
        self.assertEqual(score, 1.0)
        self.assertIn("OncoKB", reason)

        # Same variant without OncoKB call -> VUS.
        call, reason, score = classify_snv("missense_variant", gene_kb, "Y367C")
        self.assertEqual(call, "VUS")

        # On-hotspot missense always DRIVER, with or without OncoKB.
        call, reason, score = classify_snv("missense_variant", gene_kb, "V550L")
        self.assertEqual(call, "DRIVER")


class TestPhase1Integration(unittest.TestCase):
    """Black-box: phase 1 with the curated annotator must produce the same TSV
    as the v0.11 inline parsing did. A regression net for the refactor."""

    def test_toy_fgfr4_via_curated(self):
        import subprocess
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "p1.tsv"
            subprocess.run([
                "python3", "bin/phase1_annotate.py",
                "--vcf", "tests/data/toy_fgfr4_only.vcf",
                "--cna", "assets/empty.cna.tsv",
                "--fusion", "assets/empty.fusion.tsv",
                "--targets-kb", "assets/targets_kb.tsv",
                "--sample-id", "TOY_FGFR4",
                "--out", str(out),
                "--annotator", "curated",
            ], check=True, cwd=REPO_ROOT)
            content = out.read_text()
            # FGFR4 hotspot V550L must surface as DRIVER.
            self.assertIn("FGFR4", content)
            self.assertIn("V550L", content)
            self.assertIn("DRIVER", content)


if __name__ == "__main__":
    unittest.main()
