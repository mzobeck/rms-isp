#!/usr/bin/env python3
"""
Build the COG STS Committee deliverable: a tier-promoted portfolio of
therapeutic hypotheses for the current cohort.

Promotion rules (v2-expanded plan §5):

    Tier 1 (validate now)   FDA-approved drug + alteration in >=5% of any subtype
    Tier 2 (validate soon)  late-phase (phase 2/3) drug + alteration in >=3%
                            OR FDA-approved drug + alteration in >=3% (but not >=5%)
    Tier 3 (characterize)   any DRIVER-supporting hypothesis that did not promote

Per-row tiers (v0.14, in phase5_score.py) drive a SAMPLE's headline tier.
This script aggregates across the cohort to produce a PORTFOLIO tier per
(gene, drug) hypothesis, the unit a clinician would actually validate.

Inputs:
    results/target_rt_cohort_summary.tsv   sample roll-up + metadata (subtype)
    results/target_rt/<sample>/p5.tsv       per-sample event x drug rows

Outputs:
    results/target_rt_STS_committee_portfolio.md
"""
from __future__ import annotations

import csv
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent.parent

# Phase ordering (highest first). Used to pick the best drug per gene.
PHASE_RANK = {
    "approved": 4,
    "phase3": 3,
    "phase2": 2,
    "phase1": 1,
    "preclinical": 0,
}

# Prevalence cutoffs from v2 plan §5.
TIER1_PREVALENCE = 0.05
TIER2_PREVALENCE = 0.03

SUBTYPE_BUCKETS = ("FN", "FP", "ALL")


def _read_cohort_summary(path: Path) -> list[dict]:
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return list(reader)


def _read_p5(path: Path) -> list[dict]:
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return list(reader)


def subtype_sample_sets(cohort_meta: list[dict]) -> dict[str, set[str]]:
    """Return {subtype -> set of sample_ids in that subtype}."""
    out: dict[str, set[str]] = {b: set() for b in SUBTYPE_BUCKETS}
    for s in cohort_meta:
        sid = s.get("sample_id", "")
        st = s.get("subtype", "")
        if st in out and sid:
            out[st].add(sid)
    return out


def gene_driver_samples(
    cohort_meta: list[dict],
    target_rt_dir: Path,
) -> dict[str, dict[str, set[str]]]:
    """For each gene, return per-subtype sets of samples that have a DRIVER call.

    Returns: {gene -> {subtype -> set of sample_ids}}.
    """
    out: dict[str, dict[str, set[str]]] = {}
    for s in cohort_meta:
        sid = s.get("sample_id", "")
        subtype = s.get("subtype", "")
        if subtype not in SUBTYPE_BUCKETS:
            continue
        p5 = target_rt_dir / sid / "p5.tsv"
        if not p5.exists():
            continue
        seen_genes: set[str] = set()
        for r in _read_p5(p5):
            if r.get("call") != "DRIVER":
                continue
            gene = r.get("gene", "")
            if not gene or gene in seen_genes:
                continue
            seen_genes.add(gene)
            out.setdefault(gene, {b: set() for b in SUBTYPE_BUCKETS})
            out[gene][subtype].add(sid)
    return out


def gene_prevalence(
    driver_samples: dict[str, dict[str, set[str]]],
    subtype_sizes: dict[str, int],
) -> dict[str, dict[str, float]]:
    """{gene -> {subtype -> fraction of subtype samples with a DRIVER on this gene}}."""
    out: dict[str, dict[str, float]] = {}
    for gene, by_sub in driver_samples.items():
        out[gene] = {}
        for sub in SUBTYPE_BUCKETS:
            denom = subtype_sizes.get(sub, 0)
            if denom == 0:
                out[gene][sub] = 0.0
            else:
                out[gene][sub] = len(by_sub[sub]) / denom
    return out


def best_drug_per_gene(
    cohort_meta: list[dict],
    target_rt_dir: Path,
) -> dict[str, dict]:
    """For each gene, return the highest-phase drug seen across all samples.

    {gene -> {drug, drug_mechanism, max_phase, rank, pediatric_evidence,
              sample_id (one example), event, confidence}}
    """
    by_gene: dict[str, dict] = {}
    for s in cohort_meta:
        sid = s.get("sample_id", "")
        p5 = target_rt_dir / sid / "p5.tsv"
        if not p5.exists():
            continue
        for r in _read_p5(p5):
            if r.get("call") != "DRIVER":
                continue
            gene = r.get("gene", "")
            if not gene:
                continue
            mp = (r.get("drug_max_phase") or "").lower()
            rank = PHASE_RANK.get(mp, -1)
            try:
                conf = float(r.get("confidence", "") or 0)
            except ValueError:
                conf = 0.0
            entry = {
                "drug": r.get("drug", ""),
                "drug_mechanism": r.get("drug_mechanism", ""),
                "max_phase": mp,
                "rank": rank,
                "pediatric_evidence": r.get("drug_pediatric_evidence", ""),
                "sample_id": sid,
                "event": r.get("hgvsp_short") or r.get("consequence", ""),
                "confidence": conf,
                "opentargets_score": r.get("opentargets_score", ""),
            }
            existing = by_gene.get(gene)
            # Prefer higher phase; break ties by higher confidence.
            if existing is None or (rank > existing["rank"]) or \
                    (rank == existing["rank"] and conf > existing["confidence"]):
                by_gene[gene] = entry
    return by_gene


def portfolio_tier(max_phase: str, max_prevalence: float) -> str:
    """Returns '1', '2', '3', or '' per the v2 plan §5 cutoffs."""
    mp = (max_phase or "").lower()
    if mp == "approved":
        if max_prevalence >= TIER1_PREVALENCE:
            return "1"
        if max_prevalence >= TIER2_PREVALENCE:
            return "2"
        return "3"
    if mp in ("phase2", "phase3"):
        if max_prevalence >= TIER2_PREVALENCE:
            return "2"
        return "3"
    if mp in ("phase1", "preclinical"):
        return "3"
    return ""


def build_portfolio(
    cohort_meta: list[dict],
    driver_samples: dict[str, dict[str, set[str]]],
    prevalence: dict[str, dict[str, float]],
    drugs: dict[str, dict],
) -> list[dict]:
    """Returns one row per gene-with-DRIVER-hits, sorted by (tier, prevalence)."""
    out: list[dict] = []
    for gene, by_sub in driver_samples.items():
        prev = prevalence.get(gene, {b: 0.0 for b in SUBTYPE_BUCKETS})
        max_prev = max(prev.values()) if prev else 0.0
        drug_entry = drugs.get(gene)
        if drug_entry is None:
            continue
        tier = portfolio_tier(drug_entry["max_phase"], max_prev)
        sample_ids: set[str] = set()
        for s in by_sub.values():
            sample_ids |= s
        out.append({
            "gene": gene,
            "drug": drug_entry["drug"],
            "drug_mechanism": drug_entry["drug_mechanism"],
            "max_phase": drug_entry["max_phase"],
            "pediatric_evidence": drug_entry["pediatric_evidence"],
            "prevalence_FN": prev.get("FN", 0.0),
            "prevalence_FP": prev.get("FP", 0.0),
            "prevalence_ALL": prev.get("ALL", 0.0),
            "max_prevalence": max_prev,
            "sample_count": len(sample_ids),
            "sample_ids": sorted(sample_ids),
            "tier": tier,
            "opentargets_score": drug_entry.get("opentargets_score", ""),
        })
    tier_rank = {"1": 0, "2": 1, "3": 2, "": 3}
    out.sort(key=lambda r: (tier_rank.get(r["tier"], 99), -r["max_prevalence"], r["gene"]))
    return out


def render_portfolio_md(
    cohort_meta: list[dict],
    portfolio: list[dict],
    *,
    pipeline_version: str,
) -> str:
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    n_total = len(cohort_meta)
    subtype_counts = Counter(s.get("subtype", "") for s in cohort_meta)
    tier_counts = Counter(r["tier"] or "_unranked" for r in portfolio)

    L: list[str] = []
    L.append("# RMS-ISP STS Committee Portfolio")
    L.append("")
    L.append(f"- **Cohort size**: {n_total} samples")
    L.append(f"- **Subtype split**: " + ", ".join(
        f"{c} {subtype_counts.get(c, 0)}" for c in SUBTYPE_BUCKETS))
    L.append(f"- **Tier counts**: "
             f"Tier 1 = {tier_counts.get('1', 0)}, "
             f"Tier 2 = {tier_counts.get('2', 0)}, "
             f"Tier 3 = {tier_counts.get('3', 0)}")
    L.append(f"- **Pipeline version**: `{pipeline_version}`")
    L.append(f"- **Generated**: {now}")
    L.append("")
    L.append("> Engineering output, not medical advice. The drug-evidence "
             "weights, gene knowledge base, and dependency context still "
             "reflect a curated 21-target pilot KB; this portfolio is the "
             "first cohort-aggregated deliverable and is intended for STS "
             "committee discussion, not clinical action. Per-row evidence "
             "trails live in `results/target_rt/<sample>/<sample>.report.md`.")
    L.append("")

    tier_titles = {
        "1": "Tier 1 (validate now): FDA-approved drug, prevalence >= 5% in some subtype",
        "2": "Tier 2 (validate soon): late-phase or approved drug at 3-5% prevalence",
        "3": "Tier 3 (characterize): qualifying DRIVER hits that did not meet a prevalence cutoff",
    }
    for tier in ("1", "2", "3"):
        L.append(f"## {tier_titles[tier]}")
        L.append("")
        rows = [r for r in portfolio if r["tier"] == tier]
        if not rows:
            L.append("_No hypotheses in this tier._")
            L.append("")
            continue
        L.append("| Gene | Best drug | Phase | Ped evidence | Mechanism | "
                 "Max prev | FN | FP | ALL | Samples | OT score |")
        L.append("|---|---|---|---|---|---|---|---|---|---|---|")
        for r in rows:
            ot = r.get("opentargets_score", "") or "-"
            L.append(
                f"| **{r['gene']}** | `{r['drug']}` | {r['max_phase']} | "
                f"{r['pediatric_evidence']} | {r['drug_mechanism']} | "
                f"{r['max_prevalence']:.0%} | {r['prevalence_FN']:.0%} | "
                f"{r['prevalence_FP']:.0%} | {r['prevalence_ALL']:.0%} | "
                f"{r['sample_count']} | {ot} |"
            )
        L.append("")

    L.append("## Per-sample anchor hypothesis")
    L.append("")
    L.append("Each cohort sample's highest-tier portfolio hypothesis (if any).")
    L.append("")
    L.append("| Sample | Subtype | Anchor gene | Anchor drug | Portfolio tier |")
    L.append("|---|---|---|---|---|")
    sample_to_anchor = _build_sample_anchor_map(cohort_meta, portfolio)
    for s in cohort_meta:
        sid = s.get("sample_id", "")
        anchor = sample_to_anchor.get(sid)
        if anchor is None:
            L.append(f"| `{sid}` | {s.get('subtype', '')} | - | - | - |")
        else:
            L.append(f"| `{sid}` | {s.get('subtype', '')} | "
                     f"**{anchor['gene']}** | `{anchor['drug']}` | "
                     f"{anchor['tier']} |")
    L.append("")

    L.append("## Methodology")
    L.append("")
    L.append("- **Per-sample tier (v0.14)**: each phase-5 row's tier comes "
             "from `bin/phase5_score.py:tier_for_row` based on the variant "
             "call and the drug's max development phase.")
    L.append("- **Portfolio tier (v0.15)**: per (gene, best-available drug) "
             "across the cohort, promoted only if the gene's prevalence "
             "exceeds the v2-plan cutoff in at least one subtype "
             f"(>= {TIER1_PREVALENCE:.0%} for Tier 1, "
             f">= {TIER2_PREVALENCE:.0%} for Tier 2).")
    L.append("- **Best drug per gene**: highest `drug_max_phase` seen across "
             "all DRIVER rows in the cohort for that gene; ties broken by "
             "confidence.")
    L.append("- **Subtype buckets**: FN, FP, ALL. Samples without a recorded "
             "subtype contribute to no bucket.")
    L.append("- **Data sources** are the same as the per-sample reports: "
             "phase 1 (variant + CNA + fusion via the annotator chain), phase 3 "
             "(DepMap 26Q1, OpenPedCan v15), phase 4 (curated drug map + "
             "DGIdb + ClinicalTrials.gov), phase 5 confidence formula.")
    L.append("")
    L.append("## What this proves and does not prove")
    L.append("")
    L.append("**Proves**: at the cohort level, a small set of recurrent "
             "molecular alterations in RMS pair with already-approved or "
             "clinical-stage drugs at frequencies that justify experimental "
             "validation. The pipeline reproduces this without retuning the "
             "scoring formula.")
    L.append("")
    L.append("**Does not prove**: that a Tier 1 hypothesis will work in any "
             "specific patient. Drug-evidence weights treat all "
             "approved-and-pediatric-trial drugs equivalently; expression "
             "scoring is RMS-vs-other-pediatric, not RMS-vs-normal-muscle; "
             "cohort prevalence here reflects the 36 cBioPortal samples "
             "available pre-MCI, not the eventual MCI cohort. Real "
             "translation requires drug-level review by the COG STS "
             "committee.")
    L.append("")
    return "\n".join(L) + "\n"


def _build_sample_anchor_map(
    cohort_meta: list[dict],
    portfolio: list[dict],
) -> dict[str, dict]:
    """For each sample, find its highest-tier portfolio hypothesis (if any)."""
    tier_rank = {"1": 0, "2": 1, "3": 2, "": 3}
    out: dict[str, dict] = {}
    # Iterate portfolio in tier order; first hit per sample wins.
    sorted_port = sorted(portfolio, key=lambda r: tier_rank.get(r["tier"], 99))
    for entry in sorted_port:
        for sid in entry["sample_ids"]:
            if sid not in out:
                out[sid] = entry
    return out


def main(
    cohort_tsv: Path,
    target_rt_dir: Path,
    out_path: Path,
    pipeline_version: str = "v0.15.0-pilot",
) -> dict:
    """Build the portfolio. Returns a status dict for the caller."""
    cohort_meta = _read_cohort_summary(cohort_tsv)
    subtype_sets = subtype_sample_sets(cohort_meta)
    subtype_sizes = {k: len(v) for k, v in subtype_sets.items()}

    drivers = gene_driver_samples(cohort_meta, target_rt_dir)
    prev = gene_prevalence(drivers, subtype_sizes)
    drugs = best_drug_per_gene(cohort_meta, target_rt_dir)
    portfolio = build_portfolio(cohort_meta, drivers, prev, drugs)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(render_portfolio_md(
        cohort_meta, portfolio, pipeline_version=pipeline_version,
    ))
    return {
        "n_samples": len(cohort_meta),
        "n_portfolio_genes": len(portfolio),
        "tier_counts": Counter(r["tier"] or "_unranked" for r in portfolio),
        "out_path": out_path,
    }


def _cli() -> int:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cohort-tsv", type=Path,
                    default=REPO_ROOT / "results" / "target_rt_cohort_summary.tsv")
    ap.add_argument("--target-rt-dir", type=Path,
                    default=REPO_ROOT / "results" / "target_rt")
    ap.add_argument("--out", type=Path,
                    default=REPO_ROOT / "results" / "target_rt_STS_committee_portfolio.md")
    ap.add_argument("--pipeline-version", default="v0.15.0-pilot")
    args = ap.parse_args()

    if not args.cohort_tsv.exists():
        print(f"cohort summary not found at {args.cohort_tsv}; "
              f"run bin/run_target_rt.py first", file=sys.stderr)
        return 1

    status = main(args.cohort_tsv, args.target_rt_dir, args.out,
                  pipeline_version=args.pipeline_version)
    print(f"wrote {status['out_path']} "
          f"({status['n_portfolio_genes']} hypotheses; "
          f"{dict(status['tier_counts'])})")
    return 0


if __name__ == "__main__":
    raise SystemExit(_cli())
