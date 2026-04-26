"""OncoKB second-pass annotator.

Augments a prior annotation pass (typically VEP REST) with OncoKB's
oncogenic / mutation-effect classification per variant per disease.

OncoKB endpoint: POST https://www.oncokb.org/api/v1/annotate/mutations/byProteinChange
Auth:           Bearer token in Authorization header. Free for academic use:
                https://www.oncokb.org/account/register
Tumor type:     OncoTree code or generic name. RMS-relevant codes: RMS, ERMS, ARMS, PLRMS.
Cache:          One file per (gene, hgvsp_short, tumor_type) under data/oncokb_cache/.
"""
from __future__ import annotations

import json
import os
import sys
import urllib.error
import urllib.request
from pathlib import Path

from .base import (
    AugmentingAnnotator,
    Variant,
    VariantAnnotation,
    merge_annotations,
)


ONCOKB_ENDPOINT = "https://www.oncokb.org/api/v1/annotate/mutations/byProteinChange"
BATCH_SIZE = 100


def _cache_key(gene: str, alteration: str, tumor_type: str) -> str:
    safe_alt = alteration.replace("*", "X").replace("/", "_") or "none"
    return f"{gene}__{safe_alt}__{tumor_type}.json"


class OncoKBAnnotator(AugmentingAnnotator):
    """Lifts variant annotations with OncoKB's oncogenic call.

    Only annotates variants that already have a `gene` and `hgvsp_short` from
    the prior pass. Variants without those fields pass through unchanged.

    If `token` is None and ONCOKB_TOKEN is unset, the annotator logs a warning
    once on first call and degrades to a pass-through (no OncoKB enrichment).
    Cohort runs continue with whatever the prior pass produced.
    """

    def __init__(
        self,
        cache_dir: Path,
        token: str | None = None,
        tumor_type: str = "RMS",
        reference_genome: str = "GRCh38",
    ):
        self.cache_dir = Path(cache_dir)
        self.token = token or os.environ.get("ONCOKB_TOKEN", "")
        self.tumor_type = tumor_type
        self.reference_genome = reference_genome
        self._warned_no_token = False

    def annotate_with_prior(
        self,
        variants: list[Variant],
        prior: list[VariantAnnotation],
    ) -> list[VariantAnnotation]:
        if len(variants) != len(prior):
            raise ValueError("OncoKBAnnotator: variants and prior must be same length")

        if not self.token:
            if not self._warned_no_token:
                print("OncoKB: ONCOKB_TOKEN not set; passing through prior annotations",
                      file=sys.stderr)
                self._warned_no_token = True
            return list(prior)

        # Step 1: collect (idx, gene, hgvsp_short) tuples that are eligible.
        eligible: list[tuple[int, str, str]] = []
        for i, ann in enumerate(prior):
            if ann.gene and ann.hgvsp_short:
                eligible.append((i, ann.gene, ann.hgvsp_short))

        # Step 2: separate cache hits from misses.
        cache_hits: dict[int, dict] = {}
        misses: list[tuple[int, str, str]] = []
        for i, gene, alt in eligible:
            cpath = self.cache_dir / _cache_key(gene, alt, self.tumor_type)
            if cpath.exists():
                try:
                    cache_hits[i] = json.loads(cpath.read_text())
                except (json.JSONDecodeError, OSError):
                    misses.append((i, gene, alt))
            else:
                misses.append((i, gene, alt))

        # Step 3: fetch misses in batches.
        responses: dict[int, dict] = dict(cache_hits)
        for start in range(0, len(misses), BATCH_SIZE):
            chunk = misses[start:start + BATCH_SIZE]
            batch_resp = self._fetch_batch([(g, a) for _, g, a in chunk])
            for (idx, gene, alt), resp in zip(chunk, batch_resp):
                if resp is not None:
                    responses[idx] = resp
                    self._write_cache(gene, alt, resp)

        # Step 4: merge into the prior annotations.
        out: list[VariantAnnotation] = []
        for i, base in enumerate(prior):
            resp = responses.get(i)
            if resp is None:
                out.append(base)
                continue
            oncogenic = resp.get("oncogenic", "") or ""
            me = resp.get("mutationEffect") or {}
            mutation_effect = me.get("knownEffect", "") if isinstance(me, dict) else ""
            out.append(merge_annotations(
                base,
                oncogenic=oncogenic,
                mutation_effect=mutation_effect,
                source_suffix="oncokb",
            ))
        return out

    def _fetch_batch(self, queries: list[tuple[str, str]]) -> list[dict | None]:
        """POST one batch to OncoKB. Returns one response per input or None on error."""
        if not queries:
            return []
        body = json.dumps([
            {
                "id": f"q{i}",
                "hugoSymbol": gene,
                "alteration": alt,
                "tumorType": self.tumor_type,
                "referenceGenome": self.reference_genome,
            }
            for i, (gene, alt) in enumerate(queries)
        ]).encode()
        req = urllib.request.Request(
            ONCOKB_ENDPOINT,
            data=body,
            headers={
                "Content-Type": "application/json",
                "Accept": "application/json",
                "Authorization": f"Bearer {self.token}",
            },
            method="POST",
        )
        try:
            with urllib.request.urlopen(req, timeout=60) as resp:
                payload = json.loads(resp.read().decode())
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as exc:
            print(f"OncoKB batch failed ({exc}); pass-through for {len(queries)} variants",
                  file=sys.stderr)
            return [None] * len(queries)
        except json.JSONDecodeError as exc:
            print(f"OncoKB batch returned non-JSON ({exc}); pass-through", file=sys.stderr)
            return [None] * len(queries)

        if not isinstance(payload, list):
            return [None] * len(queries)

        # OncoKB returns one entry per input in input order.
        results: list[dict | None] = []
        for idx in range(len(queries)):
            if idx < len(payload) and isinstance(payload[idx], dict):
                results.append(payload[idx])
            else:
                results.append(None)
        return results

    def _write_cache(self, gene: str, alt: str, response: dict) -> None:
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        cpath = self.cache_dir / _cache_key(gene, alt, self.tumor_type)
        try:
            cpath.write_text(json.dumps(response))
        except OSError as exc:
            print(f"OncoKB cache write failed for {cpath}: {exc}", file=sys.stderr)
