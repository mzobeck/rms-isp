"""Annotator protocol + shared types + factory."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Protocol


@dataclass(frozen=True)
class Variant:
    """A single VCF row, normalized to the fields the annotators consume."""
    chrom: str
    pos: int
    ref: str
    alt: str
    info: dict[str, str] = field(default_factory=dict)


@dataclass(frozen=True)
class VariantAnnotation:
    """Per-variant annotation output. Empty consequence means 'no annotation'.

    classify_snv() in bin/phase1_annotate.py routes empty consequences to
    PASSENGER, so backends can return a degraded annotation (empty consequence,
    best-effort gene/hgvsp_short) without raising.
    """
    gene: str
    consequence: str
    hgvsp_short: str
    source: str


class Annotator(Protocol):
    """All annotators implement this single batched method."""

    def annotate_batch(self, variants: list[Variant]) -> list[VariantAnnotation]:
        ...


def get_annotator(name: str, **opts: Any) -> Annotator:
    """Factory.

    Accepted names:
        'curated'   reads gene/consequence/hgvsp from VCF INFO (toy fixtures)
        'vep_rest'  Ensembl VEP REST API (real cohorts)

    Recognized opts:
        cache_dir   Path, used by 'vep_rest'. Defaults to <repo>/data/vep_cache.
        disease     str, currently unused, reserved for future OncoKB/AlphaMissense backends.
    """
    if name == "curated":
        from .curated_vcf import CuratedVCFAnnotator
        return CuratedVCFAnnotator()
    if name == "vep_rest":
        from .vep_rest import VEPRestAnnotator
        cache_dir = opts.get("cache_dir")
        if cache_dir is None:
            cache_dir = Path(__file__).resolve().parents[2] / "data" / "vep_cache"
        return VEPRestAnnotator(cache_dir=Path(cache_dir))
    raise ValueError(f"unknown annotator: {name!r} "
                     f"(accepted: 'curated', 'vep_rest')")
