#!/usr/bin/env python3
"""
Fetch pediatric rhabdomyosarcoma trials from ClinicalTrials.gov v2 API and
build a drug -> RMS-trial-evidence cache that Phase 4 can use to upgrade the
pediatric_evidence flag on drugs that are actually being studied in pediatric
RMS patients.

Source: https://clinicaltrials.gov/api/v2/studies. Public, no auth, no rate
limit issues for our query volume.

Strategy:
    1. Query all studies with condition "rhabdomyosarcoma" and pediatric age range
    2. Walk paginated results, extract per-study (drug, phase, status, NCT)
    3. Aggregate per drug: trial_count, max_phase, any_recruiting, example_nct
    4. Write assets/ctgov_rms_drugs.tsv

Phase 4 reads this file (when present) and applies an upgrade rule:
  - drug found in pediatric RMS trial -> pediatric_evidence becomes "yes_trial"
    (unless already "yes_approved", which outranks)
  - drug NOT found in pediatric RMS trial -> no change

This makes pediatric_evidence reflect actual real-world trial activity for
RMS specifically, not just our hand-curated guess.

Usage:
    bin/fetch_clinicaltrials.py
    bin/fetch_clinicaltrials.py --condition "rhabdomyosarcoma" --max-pages 5
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent.parent
CTGOV_URL = "https://clinicaltrials.gov/api/v2/studies"

# CT.gov phase strings, ordered low -> high so we can pick the max
PHASE_ORDER = ["EARLY_PHASE1", "PHASE1", "PHASE1_PHASE2", "PHASE2", "PHASE2_PHASE3", "PHASE3", "PHASE4", "NA"]
PHASE_RANK = {p: i for i, p in enumerate(PHASE_ORDER)}

OUT_COLS = [
    "drug", "n_rms_trials", "max_phase_in_rms", "any_recruiting",
    "any_pediatric", "example_nct", "example_title",
]


def fetch_page(condition: str, page_token: str | None, page_size: int = 100) -> dict:
    params = {
        "query.cond": condition,
        "pageSize": str(page_size),
        "fields": ",".join([
            "NCTId", "BriefTitle", "Phase", "OverallStatus",
            "InterventionName", "InterventionType", "StdAge",
        ]),
    }
    if page_token:
        params["pageToken"] = page_token
    url = f"{CTGOV_URL}?{urllib.parse.urlencode(params)}"
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=60) as resp:
        return json.loads(resp.read())


def normalize_drug(name: str) -> str:
    """Match the casing convention used in our curated drug map."""
    return name.strip().lower()


# Names to skip: too generic to be useful for therapeutic-hypothesis ranking.
SKIP_DRUGS = {
    "placebo", "saline", "control", "standard of care", "best supportive care",
    "vehicle", "no intervention", "observation", "questionnaire", "blood draw",
    "biospecimen collection",
}


def extract_drug_rows(study: dict) -> list[dict]:
    """Pull per-drug evidence rows from one CT.gov study."""
    proto = study.get("protocolSection", {})
    nct = proto.get("identificationModule", {}).get("nctId", "")
    title = proto.get("identificationModule", {}).get("briefTitle", "")
    status = proto.get("statusModule", {}).get("overallStatus", "")
    phases = proto.get("designModule", {}).get("phases", []) or []
    interventions = proto.get("armsInterventionsModule", {}).get("interventions", []) or []
    ages = proto.get("eligibilityModule", {}).get("stdAges", []) or []
    is_pediatric = "CHILD" in ages
    rows: list[dict] = []
    for inter in interventions:
        if inter.get("type") not in {"DRUG", "BIOLOGICAL"}:
            continue
        name = inter.get("name") or ""
        norm = normalize_drug(name)
        if not norm or norm in SKIP_DRUGS:
            continue
        for phase in (phases or [None]):
            rows.append({
                "drug": norm,
                "phase": phase or "NA",
                "status": status,
                "nct": nct,
                "title": title,
                "is_pediatric": is_pediatric,
            })
    return rows


def aggregate(rows: list[dict]) -> list[dict]:
    """Collapse per-trial rows to per-drug summary rows."""
    by_drug: dict[str, dict] = {}
    for r in rows:
        d = r["drug"]
        agg = by_drug.setdefault(d, {
            "drug": d,
            "trials": set(),
            "phases": set(),
            "any_recruiting": False,
            "any_pediatric": False,
            "example_nct": "",
            "example_title": "",
        })
        agg["trials"].add(r["nct"])
        agg["phases"].add(r["phase"])
        if r["status"].upper() in {"RECRUITING", "ENROLLING_BY_INVITATION", "ACTIVE_NOT_RECRUITING"}:
            agg["any_recruiting"] = True
        if r["is_pediatric"]:
            agg["any_pediatric"] = True
        # Prefer the most recent or the first pediatric trial as the example.
        if not agg["example_nct"] or (r["is_pediatric"] and "pediatric" not in agg["example_title"].lower()):
            agg["example_nct"] = r["nct"]
            agg["example_title"] = r["title"]

    out: list[dict] = []
    for d, agg in by_drug.items():
        max_phase = max(agg["phases"], key=lambda p: PHASE_RANK.get(p, -1))
        out.append({
            "drug": d,
            "n_rms_trials": len(agg["trials"]),
            "max_phase_in_rms": max_phase,
            "any_recruiting": "1" if agg["any_recruiting"] else "0",
            "any_pediatric": "1" if agg["any_pediatric"] else "0",
            "example_nct": agg["example_nct"],
            "example_title": agg["example_title"][:120],
        })
    out.sort(key=lambda r: (-r["n_rms_trials"], r["drug"]))
    return out


def write_tsv(rows: list[dict], out_path: Path, *, condition: str, n_studies: int) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    with out_path.open("w", newline="") as fh:
        fh.write("## Pediatric RMS clinical trials drug map (auto-generated; do NOT hand-edit)\n")
        fh.write(f"## Source: {CTGOV_URL}\n")
        fh.write(f"## Condition query: {condition}\n")
        fh.write(f"## Generated: {timestamp}\n")
        fh.write(f"## Studies scanned: {n_studies}\n")
        fh.write(f"## Drug rows: {len(rows)}\n")
        fh.write("## Schema: drug\\tn_rms_trials\\tmax_phase_in_rms\\tany_recruiting\\tany_pediatric\\texample_nct\\texample_title\n")
        fh.write("##\n")
        fh.write("## Phase 4 reads this file (when present) and upgrades pediatric_evidence to\n")
        fh.write("## 'yes_trial' for any drug appearing here with any_pediatric=1, unless its\n")
        fh.write("## existing pediatric_evidence is already 'yes_approved' (which outranks).\n")
        fh.write("## To regenerate: bin/fetch_clinicaltrials.py\n")
        w = csv.DictWriter(fh, fieldnames=OUT_COLS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--condition", default="rhabdomyosarcoma")
    ap.add_argument("--out", type=Path, default=REPO_ROOT / "assets" / "ctgov_rms_drugs.tsv")
    ap.add_argument("--max-pages", type=int, default=20,
                    help="Cap on paginated CT.gov fetches (each page = 100 studies). 20 = 2000 studies max, plenty for RMS.")
    ap.add_argument("--page-size", type=int, default=100)
    args = ap.parse_args()

    print(f"querying CT.gov for condition='{args.condition}' (page_size={args.page_size}, max_pages={args.max_pages})", file=sys.stderr)
    all_rows: list[dict] = []
    n_studies = 0
    page_token: str | None = None
    for i in range(args.max_pages):
        page = fetch_page(args.condition, page_token, args.page_size)
        studies = page.get("studies", [])
        n_studies += len(studies)
        for s in studies:
            all_rows.extend(extract_drug_rows(s))
        page_token = page.get("nextPageToken")
        if not page_token:
            break
        print(f"  page {i+1}: {len(studies)} studies (running total {n_studies})", file=sys.stderr)
    print(f"scanned {n_studies} studies; extracted {len(all_rows)} per-trial drug rows", file=sys.stderr)

    aggregated = aggregate(all_rows)
    print(f"aggregated to {len(aggregated)} unique drugs", file=sys.stderr)

    write_tsv(aggregated, args.out, condition=args.condition, n_studies=n_studies)
    print(f"wrote {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
