"""Reads gene/consequence/hgvsp_short from VCF INFO (current-behavior backend).

Equivalent to the inline parsing that lived in bin/phase1_annotate.py before
v0.12. Toy fixtures and any cohort whose fetcher pre-bakes GENE/CONSEQUENCE
into the VCF INFO field use this backend.
"""
from __future__ import annotations

import re

from .base import Annotator, Variant, VariantAnnotation


HOTSPOT_RE = re.compile(r"(?<![A-Za-z0-9])([A-Z]\d+[A-Z*X])(?![A-Za-z0-9])")


def extract_hotspot(note: str) -> str:
    """Pull a single-letter HGVSp short (e.g., 'V550L') out of a free-text note."""
    m = HOTSPOT_RE.search(note or "")
    return m.group(1) if m else ""


class CuratedVCFAnnotator(Annotator):
    def annotate_batch(self, variants: list[Variant]) -> list[VariantAnnotation]:
        out: list[VariantAnnotation] = []
        for v in variants:
            gene = v.info.get("GENE", "")
            consequence = v.info.get("CONSEQUENCE", "")
            note = v.info.get("NOTE", "")
            out.append(VariantAnnotation(
                gene=gene,
                consequence=consequence,
                hgvsp_short=extract_hotspot(note),
                source="curated",
            ))
        return out
