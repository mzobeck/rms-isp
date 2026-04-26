"""Ensembl VEP REST API backend.

Public endpoint, no install, no 50 GB cache. Per-variant responses cache to
disk (under <repo>/data/vep_cache/) so the second run on the same cohort hits
local files only.

Endpoint: https://rest.ensembl.org/vep/human/region (POST, JSON in/out).
Rate limit: 15 req/sec, 55,000 req/hour publicly. Batches up to 200 variants
per request, so a 36-sample / 1000-sample cohort is one or two batches.
"""
from __future__ import annotations

import json
import sys
import urllib.error
import urllib.request
from pathlib import Path

from .base import Annotator, Variant, VariantAnnotation


VEP_ENDPOINT = "https://rest.ensembl.org/vep/human/region"
BATCH_SIZE = 200

# 3-letter amino acid -> 1-letter, for converting VEP's hgvsp into our short form.
AA_3TO1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Glu": "E", "Gln": "Q", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Ter": "*",   # stop codon
}


def _cache_key(v: Variant) -> str:
    return f"{v.chrom}_{v.pos}_{v.ref}_{v.alt}.json"


def hgvsp_short(hgvsp: str) -> str:
    """Convert ENST...:p.Val550Leu (or just p.Val550Leu) to V550L.

    Returns empty string for anything that does not match the simple
    single-residue substitution pattern (e.g., insertions, frameshifts).
    """
    if not hgvsp or "p." not in hgvsp:
        return ""
    p = hgvsp.split("p.", 1)[1]
    # Expect: <3letter_ref><pos><3letter_alt>, e.g. Val550Leu, Arg175His, Ter249Cys
    if len(p) < 7:
        return ""
    ref3 = p[:3]
    alt3 = p[-3:]
    middle = p[3:-3]
    if not middle.isdigit():
        return ""
    if ref3 not in AA_3TO1 or alt3 not in AA_3TO1:
        return ""
    return f"{AA_3TO1[ref3]}{middle}{AA_3TO1[alt3]}"


def pick_canonical(transcript_consequences: list[dict]) -> dict | None:
    """Return the canonical transcript_consequence object, or None."""
    if not transcript_consequences:
        return None
    for tc in transcript_consequences:
        if tc.get("canonical") == 1:
            return tc
    return transcript_consequences[0]


def parse_vep_response(obj: dict) -> tuple[str, str, str]:
    """Returns (gene, consequence, hgvsp_short) from a single VEP response object."""
    tc = pick_canonical(obj.get("transcript_consequences", []))
    if tc is None:
        # Fall back to the variant-level summary.
        return ("", obj.get("most_severe_consequence", ""), "")
    gene = tc.get("gene_symbol", "")
    cons_terms = tc.get("consequence_terms") or []
    consequence = cons_terms[0] if cons_terms else obj.get("most_severe_consequence", "")
    return (gene, consequence, hgvsp_short(tc.get("hgvsp", "")))


def _to_vep_input(v: Variant) -> str:
    """One-line VCF-style input for the VEP REST batch endpoint."""
    return f"{v.chrom} {v.pos} . {v.ref} {v.alt} . . ."


class VEPRestAnnotator(Annotator):
    def __init__(self, cache_dir: Path):
        self.cache_dir = Path(cache_dir)

    def annotate_batch(self, variants: list[Variant]) -> list[VariantAnnotation]:
        # Step 1: load cache hits, identify cache misses.
        responses: dict[int, dict] = {}
        misses: list[tuple[int, Variant]] = []
        for i, v in enumerate(variants):
            cache_path = self.cache_dir / _cache_key(v)
            if cache_path.exists():
                try:
                    responses[i] = json.loads(cache_path.read_text())
                except (json.JSONDecodeError, OSError):
                    misses.append((i, v))
            else:
                misses.append((i, v))

        # Step 2: fetch misses in batches.
        for start in range(0, len(misses), BATCH_SIZE):
            chunk = misses[start:start + BATCH_SIZE]
            chunk_responses = self._fetch_batch([v for _, v in chunk])
            for (i, v), resp in zip(chunk, chunk_responses):
                if resp is not None:
                    responses[i] = resp
                    self._write_cache(v, resp)

        # Step 3: parse + assemble in input order.
        out: list[VariantAnnotation] = []
        for i, v in enumerate(variants):
            resp = responses.get(i)
            if resp is None:
                out.append(VariantAnnotation(
                    gene="", consequence="", hgvsp_short="", source="vep_rest"))
                continue
            gene, consequence, hgvsp = parse_vep_response(resp)
            out.append(VariantAnnotation(
                gene=gene, consequence=consequence,
                hgvsp_short=hgvsp, source="vep_rest"))
        return out

    def _fetch_batch(self, variants: list[Variant]) -> list[dict | None]:
        """POST one batch to VEP REST. Returns one response per input or None on error."""
        if not variants:
            return []
        body = json.dumps({"variants": [_to_vep_input(v) for v in variants]}).encode()
        req = urllib.request.Request(
            VEP_ENDPOINT,
            data=body,
            headers={"Content-Type": "application/json",
                     "Accept": "application/json"},
            method="POST",
        )
        try:
            with urllib.request.urlopen(req, timeout=60) as resp:
                payload = json.loads(resp.read().decode())
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as exc:
            print(f"VEP REST batch failed ({exc}); degrading to empty annotations",
                  file=sys.stderr)
            return [None] * len(variants)
        except json.JSONDecodeError as exc:
            print(f"VEP REST batch returned non-JSON ({exc}); degrading",
                  file=sys.stderr)
            return [None] * len(variants)

        if not isinstance(payload, list):
            return [None] * len(variants)

        # VEP returns one entry per input variant in input order.
        # Defensive: pad/truncate to the input length.
        results: list[dict | None] = []
        for idx in range(len(variants)):
            if idx < len(payload) and isinstance(payload[idx], dict):
                results.append(payload[idx])
            else:
                results.append(None)
        return results

    def _write_cache(self, v: Variant, response: dict) -> None:
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        cache_path = self.cache_dir / _cache_key(v)
        try:
            cache_path.write_text(json.dumps(response))
        except OSError as exc:
            print(f"VEP cache write failed for {cache_path}: {exc}", file=sys.stderr)
