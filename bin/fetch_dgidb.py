#!/usr/bin/env python3
"""
Fetch drug-gene interactions from DGIdb (https://dgidb.org/api/graphql) for the
RMS-ISP target set and write them to a drug_target_map-compatible TSV.

The output TSV has the same columns as assets/drug_target_map.tsv so Phase 4
can read it directly via --drug-map-extra. Curated rows always win on
(gene, drug) collisions because their mechanism strings (e.g. CDK4/6_inhibitor,
MEK_inhibitor) are what the case-study scorecard asserts on; DGIdb only
labels interactions as generic types like inhibitor / antagonist / modulator.

Pediatric evidence cannot be inferred from DGIdb so all DGIdb-sourced rows
get pediatric_evidence = none.

Usage:
    bin/fetch_dgidb.py                                # uses targets from targets_kb.tsv
    bin/fetch_dgidb.py --out assets/dgidb_drugs.tsv   # explicit output
    bin/fetch_dgidb.py --gene CDK4 --gene FGFR4       # subset for testing

Network requirement: makes ONE GraphQL POST to dgidb.org. Run manually whenever
the curated map needs refreshing; the resulting TSV is committed to the repo
so the pipeline itself stays offline.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import urllib.request
from datetime import datetime, timezone
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent.parent
DGIDB_URL = "https://dgidb.org/api/graphql"

QUERY = """
query genes($names: [String!]!) {
  genes(names: $names) {
    nodes {
      name
      longName
      interactions {
        drug { name conceptId approved }
        interactionTypes { type }
        sources { sourceDbName }
        publications { pmid }
      }
    }
  }
}
"""


OUT_COLS = ["gene", "drug", "mechanism", "max_phase", "pediatric_evidence", "notes"]


def load_targets(targets_kb: Path) -> list[str]:
    out: list[str] = []
    with targets_kb.open() as fh:
        reader = csv.DictReader((ln for ln in fh if not ln.startswith("##")), delimiter="\t")
        for row in reader:
            out.append(row["gene"].strip())
    return out


def query_dgidb(genes: list[str], timeout: int = 30) -> dict:
    body = json.dumps({"query": QUERY, "variables": {"names": genes}}).encode()
    req = urllib.request.Request(
        DGIDB_URL, data=body,
        headers={"Content-Type": "application/json", "Accept": "application/json"},
        method="POST",
    )
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        if resp.status != 200:
            raise RuntimeError(f"DGIdb returned HTTP {resp.status}")
        return json.loads(resp.read().decode())


def normalize_drug_name(name: str) -> str:
    return name.strip().lower().replace(" ", "_")


DROPPED_TYPES_FROM_NOISE = "rows with no DGIdb interactionType (catch-all 'interaction') are dropped because they flood the report with literature-cited but non-therapeutic drug-gene pairs (e.g. doxorubicin/MYOD1)."

# Sources that are clinical-trial-listing databases rather than mechanistic
# annotations. When an interaction has only such a source, the 'interaction' is
# usually a co-occurrence in trial enrollment criteria (drug-X-given-to-patients-
# with-gene-Y-mutation) rather than a mechanistic target relationship. The
# motivating false positive: TP53 'activator' granisetron (a 5-HT3 anti-emetic)
# sourced only from ClearityFoundationClinicalTrial. This filter drops only rows
# whose ENTIRE source set is in this list; co-occurrence with any other source
# keeps the row.
LOW_TRUST_SOLE_SOURCES = {
    "ClearityFoundationClinicalTrial",
}
DROPPED_TRIAL_ONLY_FROM_NOISE = (
    f"rows whose only source is one of {sorted(LOW_TRUST_SOLE_SOURCES)} are "
    "dropped because such databases list drugs and genes from the same "
    "clinical trial without implying a mechanistic interaction."
)


def is_low_trust_sole(sources: list[str]) -> bool:
    """True iff every source in the list is in LOW_TRUST_SOLE_SOURCES."""
    if not sources:
        return False
    return all(s in LOW_TRUST_SOLE_SOURCES for s in sources)


def render_rows(payload: dict) -> tuple[list[dict], int, int]:
    """Flatten DGIdb GraphQL payload into our drug_target_map schema.

    Returns (rows, n_dropped_no_mechanism, n_dropped_low_trust).
    """
    rows: list[dict] = []
    seen: set[tuple[str, str]] = set()
    n_dropped_mech = 0
    n_dropped_low_trust = 0
    for gene_node in payload["data"]["genes"]["nodes"]:
        gene = gene_node["name"]
        for inter in gene_node.get("interactions", []):
            drug_obj = inter.get("drug") or {}
            drug_name = drug_obj.get("name")
            if not drug_name:
                continue
            types = sorted({t["type"] for t in inter.get("interactionTypes", []) if t.get("type")})
            if not types:
                # No annotated mechanism. See DROPPED_TYPES_FROM_NOISE above.
                n_dropped_mech += 1
                continue
            sources = sorted({s["sourceDbName"] for s in inter.get("sources", []) if s.get("sourceDbName")})
            if is_low_trust_sole(sources):
                # Sole source is a clinical-trial listing. See DROPPED_TRIAL_ONLY_FROM_NOISE.
                n_dropped_low_trust += 1
                continue
            drug_norm = normalize_drug_name(drug_name)
            if (gene, drug_norm) in seen:
                continue
            seen.add((gene, drug_norm))
            mechanism = "/".join(types)
            approved = bool(drug_obj.get("approved"))
            max_phase = "approved" if approved else "phase1"
            pmids = sorted({p["pmid"] for p in inter.get("publications", []) if p.get("pmid")})
            note_parts = [f"dgidb:{drug_obj.get('conceptId', '')}"]
            if sources:
                note_parts.append(f"sources={','.join(sources)}")
            if pmids[:3]:
                note_parts.append(f"pmid={','.join(str(p) for p in pmids[:3])}")
            rows.append({
                "gene": gene,
                "drug": drug_name.lower(),
                "mechanism": mechanism,
                "max_phase": max_phase,
                "pediatric_evidence": "none",
                "notes": " | ".join(note_parts),
            })
    return rows, n_dropped_mech, n_dropped_low_trust


def write_tsv(rows: list[dict], out_path: Path, *, n_genes_queried: int) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    with out_path.open("w", newline="") as fh:
        fh.write("## DGIdb-sourced drug-target interactions (auto-generated; do NOT hand-edit)\n")
        fh.write(f"## Source: {DGIDB_URL}\n")
        fh.write(f"## Generated: {timestamp}\n")
        fh.write(f"## Genes queried: {n_genes_queried}\n")
        fh.write(f"## Rows: {len(rows)}\n")
        fh.write(f"## Schema: gene\\tdrug\\tmechanism\\tmax_phase\\tpediatric_evidence\\tnotes\n")
        fh.write("## To regenerate: bin/fetch_dgidb.py\n")
        fh.write("## Phase 4 unions this with assets/drug_target_map.tsv; the curated map wins on (gene, drug) collisions because its mechanism strings drive the case-study scorecard assertions.\n")
        w = csv.DictWriter(fh, fieldnames=OUT_COLS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets-kb", type=Path, default=REPO_ROOT / "assets" / "targets_kb.tsv")
    ap.add_argument("--gene", action="append", default=[],
                    help="Restrict to specific gene(s). Default: all genes from targets_kb.tsv.")
    ap.add_argument("--out", type=Path, default=REPO_ROOT / "assets" / "dgidb_drugs.tsv")
    ap.add_argument("--dry-run", action="store_true",
                    help="Print to stdout instead of writing the file.")
    args = ap.parse_args()

    genes = args.gene or load_targets(args.targets_kb)
    print(f"querying DGIdb for {len(genes)} genes: {', '.join(genes[:5])}{'...' if len(genes) > 5 else ''}", file=sys.stderr)
    payload = query_dgidb(genes)
    rows, n_dropped_mech, n_dropped_low_trust = render_rows(payload)
    print(
        f"got {len(rows)} drug-gene interaction rows from DGIdb "
        f"(dropped {n_dropped_mech} mechanism-unannotated rows; "
        f"{DROPPED_TYPES_FROM_NOISE}) "
        f"(dropped {n_dropped_low_trust} low-trust-only rows; "
        f"{DROPPED_TRIAL_ONLY_FROM_NOISE})",
        file=sys.stderr,
    )
    if args.dry_run:
        w = csv.DictWriter(sys.stdout, fieldnames=OUT_COLS, delimiter="\t")
        w.writeheader()
        w.writerows(rows[:20])
        print(f"... ({len(rows)} total)", file=sys.stderr)
    else:
        write_tsv(rows, args.out, n_genes_queried=len(genes))
        print(f"wrote {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
