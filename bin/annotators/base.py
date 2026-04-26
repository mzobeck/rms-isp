"""Annotator protocol + shared types + factory + chain composition."""
from __future__ import annotations

from dataclasses import dataclass, field, replace
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
    """Per-variant annotation output.

    Empty `consequence` means 'no annotation'. classify_snv() in
    bin/phase1_annotate.py routes empty consequences to PASSENGER, so backends
    can return a degraded annotation without raising.

    `oncogenic` and `mutation_effect` are populated by OncoKB-class backends
    that augment a prior annotation. Empty for backends that do not produce
    such calls (curated, vep_rest).
    """
    gene: str
    consequence: str
    hgvsp_short: str
    source: str
    oncogenic: str = ""        # "Oncogenic" | "Likely Oncogenic" | "Likely Neutral" | "Inconclusive" | ""
    mutation_effect: str = ""  # "Gain-of-function" | "Loss-of-function" | "Switch-of-function" | ""


class Annotator(Protocol):
    """All annotators implement this single batched method."""

    def annotate_batch(self, variants: list[Variant]) -> list[VariantAnnotation]:
        ...


class AugmentingAnnotator(Protocol):
    """A second-pass annotator that consumes prior annotations.

    OncoKB and similar disease-aware backends need gene + hgvsp_short already
    populated (typically from a VEP-class first-pass annotator) before they
    can do their work. They implement this richer protocol; the chain wrapper
    feeds them the prior pass's output.
    """

    def annotate_with_prior(
        self,
        variants: list[Variant],
        prior: list[VariantAnnotation],
    ) -> list[VariantAnnotation]:
        ...


class _Chain:
    """Wraps a list of annotators so the chain itself looks like a single Annotator.

    First annotator runs `annotate_batch`. Each subsequent annotator either:
      - implements `annotate_with_prior` (e.g., OncoKB) -> consumes prior pass
      - implements only `annotate_batch` (e.g., a second VEP) -> overrides prior
    """

    def __init__(self, members: list[Annotator]):
        self.members = members

    def annotate_batch(self, variants: list[Variant]) -> list[VariantAnnotation]:
        if not self.members:
            return [VariantAnnotation(gene="", consequence="", hgvsp_short="",
                                       source="empty_chain") for _ in variants]
        current = self.members[0].annotate_batch(variants)
        for ann in self.members[1:]:
            if hasattr(ann, "annotate_with_prior"):
                current = ann.annotate_with_prior(variants, current)
            else:
                current = ann.annotate_batch(variants)
        return current


def get_annotator(name: str, **opts: Any) -> Annotator:
    """Factory.

    Accepts a single backend name or a comma-separated chain.

    Single names:
        'curated'   reads gene/consequence/hgvsp from VCF INFO (toy fixtures)
        'vep_rest'  Ensembl VEP REST API (real cohorts)
        'oncokb'    OncoKB augmentation (requires prior gene + hgvsp_short)

    Chain example:
        'vep_rest,oncokb'   VEP first-pass, OncoKB second-pass

    Recognized opts:
        cache_dir         Path, used by 'vep_rest'. Defaults to <repo>/data/vep_cache.
        oncokb_cache_dir  Path, used by 'oncokb'. Defaults to <repo>/data/oncokb_cache.
        oncokb_token      str, OncoKB API token. Defaults to env var ONCOKB_TOKEN.
        oncokb_tumor_type str, OncoTree code or disease name. Defaults to 'RMS'.
        disease           str, generic pass-through; future backends may consume.
    """
    if "," in name:
        members = [get_annotator(part.strip(), **opts) for part in name.split(",") if part.strip()]
        return _Chain(members)

    if name == "curated":
        from .curated_vcf import CuratedVCFAnnotator
        return CuratedVCFAnnotator()
    if name == "vep_rest":
        from .vep_rest import VEPRestAnnotator
        cache_dir = opts.get("cache_dir")
        if cache_dir is None:
            cache_dir = Path(__file__).resolve().parents[2] / "data" / "vep_cache"
        return VEPRestAnnotator(cache_dir=Path(cache_dir))
    if name == "oncokb":
        from .oncokb import OncoKBAnnotator
        cache_dir = opts.get("oncokb_cache_dir")
        if cache_dir is None:
            cache_dir = Path(__file__).resolve().parents[2] / "data" / "oncokb_cache"
        token = opts.get("oncokb_token")
        tumor_type = opts.get("oncokb_tumor_type", opts.get("disease", "RMS"))
        return OncoKBAnnotator(
            cache_dir=Path(cache_dir),
            token=token,
            tumor_type=tumor_type,
        )

    raise ValueError(f"unknown annotator: {name!r} "
                     f"(accepted: 'curated', 'vep_rest', 'oncokb', "
                     f"or comma-separated chain)")


def merge_annotations(
    base: VariantAnnotation,
    *,
    oncogenic: str = "",
    mutation_effect: str = "",
    source_suffix: str = "",
) -> VariantAnnotation:
    """Helper for second-pass annotators: returns base with selected fields overlaid."""
    new_source = f"{base.source}+{source_suffix}" if source_suffix else base.source
    return replace(
        base,
        oncogenic=oncogenic or base.oncogenic,
        mutation_effect=mutation_effect or base.mutation_effect,
        source=new_source,
    )
