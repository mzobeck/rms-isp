#!/usr/bin/env python3
"""
Cohort-level visualizations for the RMS-ISP TARGET-RT runner.

Reads results/target_rt_cohort_summary.tsv plus per-sample p5.tsv outputs,
writes a long-format aggregation TSV and three SVG charts:

  - cohort_mechanisms.svg     horizontal bar chart of top_mechanism counts
  - cohort_druggability.svg   gene x subtype fraction matrix (always 21 x 4)
  - cohort_per_sample.svg     gene x sample heatmap (only when N <= 100)

Pure stdlib. Called from bin/run_target_rt.py at the end of a cohort run, or
runnable standalone for fast viz iteration without re-running 36 pipelines.
"""
from __future__ import annotations

import csv
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent

# Confidence threshold: cells below this count as "not druggable" and render blank.
# Pinned to 0.10 to match the passenger-sanity assertion in tests/cases.toml.
DRUGGABILITY_THRESHOLD = 0.10

# Color choices (committed here so the test can assert on hex values).
HEATMAP_RAMP = ["#eff3ff", "#bdd7e7", "#6baed6", "#3182bd", "#08519c"]
HEATMAP_BINS = [(0.10, 0.20), (0.20, 0.40), (0.40, 0.60), (0.60, 0.80), (0.80, 1.01)]
EMPTY_CELL = "#f5f5f5"
SUBTYPE_PALETTE = {"FN": "#377eb8", "FP": "#e41a1c", "ALL": "#984ea3"}

# Suppress per-sample heatmap above this cohort size (illegible).
PER_SAMPLE_MAX_N = 100
