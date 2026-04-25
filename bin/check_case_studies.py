#!/usr/bin/env python3
"""
RMS-ISP case-study scorecard.

Reads tests/cases.toml, runs the pipeline against each declared toy fixture,
checks each case's assertions, and writes a markdown scorecard. Exits with
status 1 if any case fails. Intended to be invoked from the repo root.

Run order per case:
    bin/phase1_annotate.py -> bin/phase2_structure.py -> bin/phase3_dependency.py
    -> bin/phase4_drugs.py -> bin/phase5_score.py

Each case's scored TSV (Phase 5 output) is loaded and the assertions are
evaluated against it. The pipeline is invoked directly (not via Nextflow) so
this script runs cleanly from a thin CI environment with only Python 3.11+.

Usage:
    bin/check_case_studies.py                     # run all cases, write scorecard
    bin/check_case_studies.py --cases case3_fgfr4 # subset
    bin/check_case_studies.py --quiet             # suppress per-step output

Outputs:
    results/scorecard.md  (human-readable)
    results/scorecard.json (machine-readable, for CI to parse)
"""
from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
import tomllib
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parent.parent


def run(cmd: list[str], *, quiet: bool) -> None:
    if not quiet:
        print(f"  $ {' '.join(cmd)}", file=sys.stderr)
    subprocess.run(cmd, check=True, cwd=REPO_ROOT)


def run_pipeline(case: dict, *, quiet: bool) -> Path:
    """Execute Phases 1-5 for one case; return the scored TSV path."""
    sid = case["sample_id"]
    workdir = REPO_ROOT / "results" / "scorecard" / sid
    workdir.mkdir(parents=True, exist_ok=True)

    vcf = case["vcf"]
    cna = case.get("cna", "assets/empty.cna.tsv")
    fusion = case.get("fusion", "assets/empty.fusion.tsv")

    p1 = workdir / "p1.tsv"
    p2 = workdir / "p2.tsv"
    p3 = workdir / "p3.tsv"
    p4 = workdir / "p4.tsv"
    p5_tsv = workdir / "p5.tsv"
    p5_md = workdir / "report.md"

    run(["python3", "bin/phase1_annotate.py",
         "--vcf", vcf, "--cna", cna, "--fusion", fusion,
         "--targets-kb", "assets/targets_kb.tsv",
         "--sample-id", sid, "--out", str(p1)], quiet=quiet)
    run(["python3", "bin/phase2_structure.py",
         "--in", str(p1), "--out", str(p2)], quiet=quiet)
    run(["python3", "bin/phase3_dependency.py",
         "--in", str(p2), "--depmap", "assets/depmap_rms_summary.tsv",
         "--subtype", case["subtype"], "--out", str(p3)], quiet=quiet)
    phase4_cmd = [
        "python3", "bin/phase4_drugs.py",
        "--in", str(p3), "--drug-map", "assets/drug_target_map.tsv",
        "--out", str(p4),
    ]
    extra = REPO_ROOT / "assets" / "dgidb_drugs.tsv"
    if extra.exists() and extra.stat().st_size > 0:
        phase4_cmd += ["--drug-map-extra", str(extra)]
    ctgov = REPO_ROOT / "assets" / "ctgov_rms_drugs.tsv"
    if ctgov.exists() and ctgov.stat().st_size > 0:
        phase4_cmd += ["--ctgov-trials", str(ctgov)]
    run(phase4_cmd, quiet=quiet)
    run(["python3", "bin/phase5_score.py",
         "--in", str(p4), "--vcf", vcf,
         "--sample-id", sid, "--subtype", case["subtype"],
         "--out-tsv", str(p5_tsv), "--out-md", str(p5_md)], quiet=quiet)
    return p5_tsv


def load_scored(path: Path) -> list[dict]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def assert_mechanism_in_top_n(rows: list[dict], mechanism: str, top_n: int) -> tuple[bool, str]:
    top = rows[:top_n]
    for r in top:
        if mechanism in (r.get("drug_mechanism") or ""):
            return True, f"found `{r['drug']}` ({r['drug_mechanism']}) at top-{top_n}"
    seen = ", ".join(f"{r['drug']}({r['drug_mechanism']})" for r in top if r.get("drug"))
    return False, f"no row in top-{top_n} matched mechanism substring `{mechanism}`. Saw: {seen or '(no drug rows)'}"


def assert_gene_event_in_top_n(rows: list[dict], gene: str, top_n: int) -> tuple[bool, str]:
    top = rows[:top_n]
    if any(r.get("gene") == gene for r in top):
        return True, f"`{gene}` present in top-{top_n}"
    seen = ", ".join(r.get("gene", "") for r in top)
    return False, f"`{gene}` not in top-{top_n}. Saw: {seen}"


def assert_passenger_below(rows: list[dict], max_conf: float) -> tuple[bool, str]:
    seen_passengers = [r for r in rows if r.get("call") == "PASSENGER"]
    if not seen_passengers:
        return False, "no PASSENGER rows in scored output"
    offenders = [r for r in seen_passengers if float(r.get("confidence", "0")) >= max_conf]
    if offenders:
        msg = "; ".join(f"{r['gene']}({r['confidence']})" for r in offenders)
        return False, f"{len(offenders)} passenger(s) at or above {max_conf}: {msg}"
    seen = ", ".join(f"{r['gene']}({r['confidence']})" for r in seen_passengers[:5])
    return True, f"all {len(seen_passengers)} passenger(s) below {max_conf}. e.g.: {seen}"


def evaluate_assertion(rows: list[dict], a: dict) -> tuple[bool, str]:
    kind = a["kind"]
    if kind == "mechanism_in_top_n":
        return assert_mechanism_in_top_n(rows, a["mechanism"], int(a["top_n"]))
    if kind == "gene_event_in_top_n":
        return assert_gene_event_in_top_n(rows, a["gene"], int(a["top_n"]))
    if kind == "passenger_below":
        return assert_passenger_below(rows, float(a["max_confidence"]))
    return False, f"unknown assertion kind `{kind}`"


def render_markdown(report: dict) -> str:
    L: list[str] = []
    overall = "PASS" if report["pass"] else "FAIL"
    L.append("# RMS-ISP Pilot Case-Study Scorecard")
    L.append("")
    L.append(f"- **Overall**: **{overall}**  ({report['n_pass']} pass, {report['n_fail']} fail of {report['n_assertions']} assertions across {report['n_cases']} cases)")
    L.append(f"- **Run timestamp (UTC)**: {report['timestamp']}")
    L.append(f"- **Pipeline version**: `{report['pipeline_version']}`")
    L.append("")
    L.append("## Per-case results")
    L.append("")
    L.append("| Case | Pilot # | Sample | Status | Assertions | Notes |")
    L.append("|---|---|---|---|---|---|")
    for c in report["cases"]:
        status = "PASS" if c["pass"] else "**FAIL**"
        notes = "; ".join(a["detail"] for a in c["assertions"] if not a["pass"]) or "ok"
        L.append(f"| `{c['id']}` | {c['case_study']} | {c['sample_id']} | {status} | {c['n_pass']}/{c['n_total']} | {notes} |")
    L.append("")
    L.append("## Detailed per-assertion results")
    L.append("")
    for c in report["cases"]:
        L.append(f"### `{c['id']}` — {c['description']}")
        L.append("")
        for a in c["assertions"]:
            mark = "PASS" if a["pass"] else "**FAIL**"
            L.append(f"- {mark} — `{a['kind']}` — {a['detail']}")
        L.append("")
    L.append("---")
    L.append("")
    L.append("Source: `tests/cases.toml`. Pipeline: see `plans/pilot_ccdi_mci/rms_translational_pilot_project.md` §3 for case definitions.")
    return "\n".join(L) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", nargs="*", help="If set, only run cases with these ids.")
    ap.add_argument("--quiet", action="store_true")
    ap.add_argument("--out-md", type=Path, default=REPO_ROOT / "results" / "scorecard.md")
    ap.add_argument("--out-json", type=Path, default=REPO_ROOT / "results" / "scorecard.json")
    args = ap.parse_args()

    cases_path = REPO_ROOT / "tests" / "cases.toml"
    with cases_path.open("rb") as fh:
        cases = tomllib.load(fh)["case"]
    if args.cases:
        cases = [c for c in cases if c["id"] in set(args.cases)]
    if not cases:
        print("no cases selected", file=sys.stderr)
        return 1

    report: dict[str, Any] = {
        "timestamp": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
        "pipeline_version": "v0.7.0-pilot",
        "n_cases": len(cases),
        "cases": [],
        "n_assertions": 0,
        "n_pass": 0,
        "n_fail": 0,
    }

    for case in cases:
        if not args.quiet:
            print(f">>> {case['id']} ({case['sample_id']}) — {case['description']}", file=sys.stderr)
        scored_tsv = run_pipeline(case, quiet=args.quiet)
        rows = load_scored(scored_tsv)
        case_results = []
        n_pass = 0
        for a in case["asserts"]:
            ok, detail = evaluate_assertion(rows, a)
            n_pass += 1 if ok else 0
            case_results.append({**a, "pass": ok, "detail": detail})
        case_pass = n_pass == len(case["asserts"])
        report["cases"].append({
            "id": case["id"],
            "case_study": case["case_study"],
            "description": case["description"],
            "sample_id": case["sample_id"],
            "pass": case_pass,
            "n_pass": n_pass,
            "n_total": len(case["asserts"]),
            "assertions": case_results,
        })
        report["n_assertions"] += len(case["asserts"])
        report["n_pass"] += n_pass
        report["n_fail"] += len(case["asserts"]) - n_pass

    report["pass"] = report["n_fail"] == 0

    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(report, indent=2) + "\n")
    args.out_md.write_text(render_markdown(report))
    if not args.quiet:
        print("", file=sys.stderr)
        print(args.out_md.read_text())

    return 0 if report["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
