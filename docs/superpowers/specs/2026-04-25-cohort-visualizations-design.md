# Cohort Visualizations Design

- **Status**: Draft, awaiting user approval
- **Author**: Claude (session 2026-04-25)
- **Targets**: v0.11.0-pilot
- **Roadmap item**: walkthrough §8 #2 ("Cohort-level visualizations")

## 1. Goal

Add static visualizations to `results/target_rt_cohort_summary.md` so readers can see cohort-level structure at a glance (mechanism mix, per-target druggability) instead of scanning a wide table. Pictures must render on GitHub directly from the markdown file. The pilot has 36 samples; the design must keep working at hundreds or thousands of samples without code changes.

## 2. Non-goals

- Interactive charts, JavaScript, or browser-only interactivity. (Standalone HTML was considered and rejected; embedding in the existing markdown wins because the report already exists and is the place readers go.)
- PNG output. (Stdlib SVG meets the same need without adding a dependency.)
- New scoring components, new gates, weight changes, or any modification to the scorecard. Visualization is downstream reporting, not part of the scientific gate.
- Mutation-level oncoprint (DRIVER / VUS / PASSENGER per gene per sample). Possible later; not in this iteration.
- Cross-cohort comparison (TARGET-RT vs MSK-IMPACT vs hypothetical future MCI). Each cohort gets its own report; cross-cohort viz is a separate spec.

## 3. Constraints

- Pure stdlib Python. The "pipeline is stdlib only" rule covers the cohort runner; viz inherits the same constraint. SVG is generated with `xml.etree.ElementTree`.
- No em dashes in committed markdown (CI grep gate).
- Must scale: design must produce a useful, legible report at N = 36, N = 200, and N = 5000 samples without code changes.
- Must be deterministic so tests can assert on SVG content.

## 4. Architecture

```
results/target_rt/<sample>/p5.tsv  (one per sample, written by run_target_rt.py)
                  |
                  v
   bin/cohort_visualize.py
                  |
   +--------------+--------------+----------------------------+
   v                             v                            v
results/target_rt_cohort_         cohort_mechanisms.svg     cohort_druggability.svg
gene_matrix.tsv                   (mechanism bar chart)     (gene x subtype heatmap)
                                                                  +
                                                            cohort_per_sample.svg
                                                            (only if N <= 100)
```

Three artifacts plus a machine-readable aggregation TSV. Each is a single SVG file under `results/target_rt/`. The cohort runner embeds them in the markdown via `<img>` tags with relative paths.

## 5. Components

### 5.1 `bin/cohort_visualize.py` (new)

CLI: `python3 bin/cohort_visualize.py [--cohort-tsv PATH] [--target-rt-dir PATH] [--out-dir PATH]`. Defaults match the canonical cohort runner output paths.

Steps:

1. **Load**. Read the cohort summary TSV (per-sample metadata, top-1 hit). For each `sample_id` in that TSV, read `<target-rt-dir>/<sample_id>/p5.tsv` (full event x drug rows). The cohort summary is the authoritative sample list, not directory listing, so any sample that produced no scored rows is correctly excluded.
2. **Aggregate**. For each (sample, gene), compute `max_confidence` across all drug rows. Write long-format `results/target_rt_cohort_gene_matrix.tsv` with columns `sample_id, study, subtype, gene, max_confidence`. Drop cells below the threshold (see §6).
3. **Render mechanism chart**. Horizontal bar chart of `top_mechanism` counts from the cohort summary. Bin everything below 1% of cohort size into a single "other (N items)" bar.
4. **Render druggability chart**. Gene rows (the 21 targets in `assets/targets_kb.tsv`, sorted by overall fraction descending) x subtype columns (`FN`, `FP`, `ALL`, `whole cohort`). Cell color = fraction of subtype samples with that gene's `max_confidence` above the threshold. Always 21 x 4 cells regardless of N. Cell text = the fraction (e.g., "0.42") since there are only 84 cells.
5. **Render per-sample heatmap conditionally**. If `N <= 100`, render a gene x sample heatmap with cells colored by `max_confidence` (0 to 1, blank below threshold). Subtype color stripe across the top. Sample text labels only when `N <= 50`. If `N > 100`, skip the file and emit a stub note for the cohort runner.

Idempotent: rerunning overwrites the SVGs.

### 5.2 `bin/run_target_rt.py` (modified)

Two changes:

1. **Order**: write the cohort TSV first; then call `cohort_visualize.main()` (in-process import, not subprocess), which returns a status dict like `{"mechanisms": Path, "druggability": Path, "per_sample": Path | None, "n_samples": int}`. Then call `write_cohort_md(rows, viz=status)`. If visualize raises, log the traceback and pass `viz=None` to `write_cohort_md` so the report still writes (the cohort TSV is the authoritative artifact and must not be blocked by a viz bug).
2. `write_cohort_md()` accepts a new optional `viz` parameter and adds three new sections after "Mechanism prevalence across the cohort":
   - "Mechanism distribution": `<img src="target_rt/cohort_mechanisms.svg" alt="Top mechanism counts across the cohort">`.
   - "Per-target cohort druggability": `<img src="target_rt/cohort_druggability.svg" alt="Per-target druggability fraction by subtype">`.
   - "Per-sample heatmap": if `viz["per_sample"]` is a Path, embed `<img src="target_rt/cohort_per_sample.svg" alt="...">`. Otherwise the section text reads "Suppressed at N=<n>; see per-target chart above for cohort-level patterns." (No em dash.)
   - If `viz` is None (visualize failed), all three sections are omitted and a single line is added: "Visualizations failed to render this run; see job log."

### 5.3 Tests: `tests/test_cohort_visualize.py` (new)

Stdlib `unittest`. Fixtures live under `tests/data/cohort_viz/`: a small synthetic cohort (3 samples, 2 genes) and a large-N cohort (150 samples, generated programmatically in test setup so it does not bloat the repo). Assertions:

- **Aggregation**: given the small fixture, the output TSV has exactly the expected sample x gene rows in the expected sort order.
- **Mechanism SVG**: well-formed XML; `<rect>` bar count equals number of distinct mechanisms in the small fixture; the "other (N items)" bar appears when the long-tail-binning rule fires.
- **Druggability SVG**: 21 x 4 = 84 cell `<rect>` elements when invoked against the full real target list (`assets/targets_kb.tsv`) plus the small fixture; cell-color hex matches the documented ramp from §6 for known cell values (e.g., a cell with fraction `0.42` gets the third bin's hex `#6baed6`).
- **Per-sample SVG threshold**: renders an SVG file for the small (N=3) fixture; returns None and writes no file for the large-N (N=150) fixture.
- **Determinism**: run the full pipeline twice on the small fixture in two temp dirs; the resulting SVG and TSV files are byte-identical.

CI: existing `case-study-scorecard` workflow runs `python3 -m unittest tests.test_cohort_visualize` before the scorecard step. No new workflow.

## 6. Threshold and color choices

- **Druggability threshold**: a (sample, gene) cell counts as "druggable" if `max_confidence >= 0.10`. This matches the existing passenger-sanity-check threshold from `tests/cases.toml`. Documented inline at the top of `bin/cohort_visualize.py` so reviewers can see and change it.
- **Color ramp** (heatmap cells): ColorBrewer Blues 5-class as fixed hex codes: `#eff3ff #bdd7e7 #6baed6 #3182bd #08519c`. Cell value bins: `[0.10,0.20)`, `[0.20,0.40)`, `[0.40,0.60)`, `[0.60,0.80)`, `[0.80,1.00]`. Picked for accessibility (single-hue, colorblind-tolerant) and because it's stdlib-renderable as fixed hex.
- **Subtype stripe** (categorical): `FN = #377eb8`, `FP = #e41a1c`, `ALL = #984ea3` (ColorBrewer Set1, three distinct colors). Documented at the top of `bin/cohort_visualize.py` as a single `PALETTE` dict.
- **Empty cells** (below threshold): rendered as `#f5f5f5` (neutral light gray) with no text, so absence is visually distinct from "low but nonzero."

## 7. Scaling behavior

| Cohort size N | Mechanism chart | Druggability chart | Per-sample heatmap |
|---|---|---|---|
| 1 to 50  | Full  | Full | Full, with sample-ID labels |
| 51 to 100 | Full | Full | Full, no sample-ID labels (subtype stripe only) |
| 101 to ~5000 | Full | Full | Suppressed; markdown shows a stub note |
| 5000+ | Full (long-tail binned) | Full | Suppressed |

The aggregation TSV is always emitted regardless of N. Anyone wanting a different visualization at any cohort size can read it directly.

## 8. Risks

- **Color choices look ugly or fail accessibility**. Mitigation: emit the SVGs locally and eyeball before committing the implementation; pick a single-hue ramp (Blues from ColorBrewer is widely cited and stdlib-renderable as fixed hex codes). If accessibility becomes a real concern, swap palette in one place.
- **Stdlib SVG generation produces non-deterministic attribute ordering across Python versions**. Mitigation: `ElementTree` writes attributes in insertion order in CPython 3.8+ which is what the project already targets; tests pin determinism on the developer's machine and CI runner. If a future Python release changes this, the test catches it.
- **Walkthrough drift**. Walkthrough §8 #2 reads "mechanism distribution and a gene-by-sample heatmap." We are delivering the mechanism distribution, replacing the per-sample heatmap with a per-target cohort druggability chart for scale, and keeping the per-sample heatmap as a conditional bonus. Walkthrough text needs a one-line update to match.

## 9. Test plan (acceptance)

- `python3 -m unittest tests.test_cohort_visualize` exits 0.
- `python3 bin/cohort_visualize.py` produces three SVG files plus the aggregation TSV in `results/target_rt/` (after `run_target_rt.py` has been run).
- `open results/target_rt_cohort_summary.md` (or view on GitHub) shows the two embedded charts inline with the existing tables.
- `python3 bin/check_case_studies.py` still exits 0 (scorecard untouched).
- The em-dash gate in `.github/workflows/ci.yml` (markdown-lint job) reports zero hits across tracked markdown.
- Synthetic large-cohort test (N=150 in fixture): per-sample SVG is omitted, markdown stub note is present.

## 10. Out of scope, captured for later

- Drug-level breakouts (top mechanism per subtype as a stacked bar, etc.).
- Mutation-call oncoprint (DRIVER / VUS / PASSENGER per gene per sample) from phase 1 output.
- Cross-cohort comparison.
- Calibrating the 0.10 cell threshold (tied to the broader confidence-calibration work, which is a v1.0 concern per walkthrough §5).
- Per-cohort viz parameterization (today the design hardcodes paths to TARGET-RT; next time another cohort is added, generalize).
