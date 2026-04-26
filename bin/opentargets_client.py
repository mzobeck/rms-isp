#!/usr/bin/env python3
"""
Thin OpenTargets v4 GraphQL client used by phase 3.

Public endpoint, no auth, no rate limits at our volume:
    https://api.platform.opentargets.org/api/v4/graphql

For a given disease (EFO ID), fetches the full table of OpenTargets-associated
targets and caches the resulting {gene_symbol -> association_score} map in a
single JSON file under data/opentargets_cache/<efo>.json. Subsequent runs hit
the cache. Per-gene lookups read the cached map.

This is a one-shot disease-side query rather than per-gene lookups: a typical
cancer EFO has 2000-3000 associated targets total, well within a single
paginated GraphQL response. Querying per-gene with the target side does not
work because the default associatedDiseases page (25 rows) often does not
include the disease of interest for genes whose top scores are with other
diseases (e.g. CDK4's top diseases are melanoma and breast cancer, not RMS).
"""
from __future__ import annotations

import json
import sys
import urllib.error
import urllib.request
from pathlib import Path

OT_ENDPOINT = "https://api.platform.opentargets.org/api/v4/graphql"

REPO_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_CACHE_DIR = REPO_ROOT / "data" / "opentargets_cache"
DEFAULT_EFO_ID = "EFO_0002918"   # rhabdomyosarcoma

# OpenTargets caps page size at 3000. RMS has ~2400 associated targets which
# fits in one page. For diseases that exceed 3000 (uncommon for cancer EFOs),
# any genes outside the top 3000 will read as score=0 in the cohort runner;
# we accept that for v0.17 and revisit with a multi-page fetcher if it bites.
DEFAULT_PAGE_SIZE = 3000

DISEASE_TARGETS_QUERY = """
query DiseaseTargets($efo: String!, $size: Int!) {
  disease(efoId: $efo) {
    id
    name
    associatedTargets(page: {index: 0, size: $size}) {
      count
      rows {
        score
        target { id approvedSymbol }
      }
    }
  }
}
"""


def _post(query: str, variables: dict, *, timeout: int = 60) -> dict | None:
    body = json.dumps({"query": query, "variables": variables}).encode()
    req = urllib.request.Request(
        OT_ENDPOINT,
        data=body,
        headers={
            "Content-Type": "application/json",
            "Accept": "application/json",
        },
        method="POST",
    )
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return json.loads(resp.read().decode())
    except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as exc:
        print(f"OpenTargets POST failed ({exc})", file=sys.stderr)
        return None
    except json.JSONDecodeError as exc:
        print(f"OpenTargets returned non-JSON ({exc})", file=sys.stderr)
        return None


def _disease_cache_path(cache_dir: Path, efo_id: str) -> Path:
    return cache_dir / f"{efo_id}.json"


def fetch_disease_target_scores(
    efo_id: str,
    cache_dir: Path,
    *,
    page_size: int = DEFAULT_PAGE_SIZE,
) -> dict | None:
    """Returns {disease_name: str, gene_scores: {symbol: score}, count: int}.

    Cached at <cache_dir>/<efo>.json. Returns None on lookup failure.
    """
    cache = _disease_cache_path(cache_dir, efo_id)
    if cache.exists():
        try:
            return json.loads(cache.read_text())
        except (json.JSONDecodeError, OSError):
            pass

    payload = _post(DISEASE_TARGETS_QUERY, {"efo": efo_id, "size": page_size})
    if payload is None:
        return None
    try:
        disease = payload["data"]["disease"]
        if disease is None:
            return None
        associated = disease.get("associatedTargets") or {}
        rows = associated.get("rows") or []
        gene_scores: dict[str, float] = {}
        for row in rows:
            target = row.get("target") or {}
            symbol = target.get("approvedSymbol", "")
            if not symbol:
                continue
            try:
                score = float(row.get("score") or 0.0)
            except (TypeError, ValueError):
                score = 0.0
            gene_scores[symbol] = score
        summary = {
            "disease_name": disease.get("name", ""),
            "efo_id": efo_id,
            "count": int(associated.get("count") or 0),
            "gene_scores": gene_scores,
        }
    except (KeyError, TypeError, ValueError) as exc:
        print(f"OpenTargets parse error for {efo_id}: {exc}", file=sys.stderr)
        return None

    cache.parent.mkdir(parents=True, exist_ok=True)
    try:
        cache.write_text(json.dumps(summary))
    except OSError as exc:
        print(f"OpenTargets cache write failed for {cache}: {exc}", file=sys.stderr)
    return summary


def lookup_gene_disease(
    gene_symbol: str,
    efo_id: str = DEFAULT_EFO_ID,
    cache_dir: Path | None = None,
) -> dict | None:
    """Returns {association_score, matched_disease_name} for a gene-disease pair.

    Returns None if the disease lookup itself fails. Returns a zero-score
    summary when the disease lookup succeeds but the gene is not among its
    associated targets (signal: OT does not list this gene for this disease).
    """
    cdir = cache_dir or DEFAULT_CACHE_DIR
    summary = fetch_disease_target_scores(efo_id, cdir)
    if summary is None:
        return None
    score = summary["gene_scores"].get(gene_symbol, 0.0)
    return {
        "association_score": score,
        "matched_disease_name": summary["disease_name"] if score > 0 else "",
        "associated_disease_count": summary["count"],
    }
