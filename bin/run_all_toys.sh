#!/usr/bin/env bash
# Run the RMS-ISP pipeline on all bundled toy patients, then generate a
# cross-sample summary. Use during development to eyeball-check that the
# pilot case studies (1=MTAP/PRMT5, 2=CDK4/CDK4-6i, 3=FGFR4/FGFRi,
# 4=RAS/MEKi) all land their expected drug classes near the top.
#
# Usage:
#     bin/run_all_toys.sh              # run from the repo root
#
# Outputs land in results/<sample_id>/ alongside a top-level
# results/multisample_summary.md.

set -euo pipefail

cd "$(dirname "$0")/.."
PROJECT_ROOT="$PWD"

declare -a SAMPLES=(
    "TOY_TUMOR:tests/data/toy_patient.vcf:assets/empty.cna.tsv:assets/empty.fusion.tsv:FN"
    "TOY_FP_CDK4:tests/data/toy_fp_cdk4amp.vcf:tests/data/toy_fp_cdk4amp.cna.tsv:tests/data/toy_fp_cdk4amp.fusion.tsv:FP"
    "TOY_MTAP_NULL:tests/data/toy_mtap_null.vcf:tests/data/toy_mtap_null.cna.tsv:tests/data/toy_mtap_null.fusion.tsv:FN"
)

mkdir -p results
SUMMARY="results/multisample_summary.md"

{
    echo "# RMS-ISP Multi-Sample Run Summary"
    echo
    echo "Run timestamp (UTC): $(date -u +'%Y-%m-%d %H:%M UTC')"
    echo
    echo "| Sample | Subtype | Top drug | Mechanism | Confidence |"
    echo "|---|---|---|---|---|"
} > "$SUMMARY"

for entry in "${SAMPLES[@]}"; do
    IFS=":" read -r sid vcf cna fusion subtype <<< "$entry"
    echo ">>> Running $sid (subtype=$subtype)"
    rm -rf "results/$sid" "work" ".nextflow" .nextflow.log* 2>/dev/null || true

    nextflow run main.nf -profile laptop \
        --input  "$vcf"  \
        --cna    "$cna"  \
        --fusion "$fusion" \
        --sample_id "$sid" \
        --subtype "$subtype" \
        --outdir "results/$sid" \
        > "results/$sid.nf.log" 2>&1

    # Pull rank-1 row from the markdown report (the row immediately under the
    # "| Rank | Gene | Event | ..." header in the top-ranked table).
    report="results/$sid/phase5/${sid}.report.md"
    top_line=$(awk '/^\| 1 \|/ {print; exit}' "$report" || true)
    if [[ -n "$top_line" ]]; then
        gene=$(echo "$top_line"      | awk -F'|' '{gsub(/[* ]/,"",$3); print $3}')
        event=$(echo "$top_line"     | awk -F'|' '{gsub(/^[ ]*|[ ]*$/,"",$4); print $4}')
        drug=$(echo "$top_line"      | awk -F'|' '{gsub(/[ `]/,"",$7); print $7}')
        mechanism=$(echo "$top_line" | awk -F'|' '{gsub(/^[ ]*|[ ]*$/,"",$8); print $8}')
        confidence=$(echo "$top_line"| awk -F'|' '{gsub(/[ *]/,"",$11); print $11}')
        echo "| $sid | $subtype | \`$drug\` ($gene $event) | $mechanism | $confidence |" >> "$SUMMARY"
    else
        echo "| $sid | $subtype | (no rank-1 found) | | |" >> "$SUMMARY"
    fi
done

{
    echo
    echo "## Pilot case-study expectations"
    echo
    echo "Per \`plans/pilot_ccdi_mci/rms_translational_pilot_project.md\` §3 case studies:"
    echo
    echo "| Case | Tumor profile | Expected top drug class | Test fixture | Eyeball-verify |"
    echo "|---|---|---|---|---|"
    echo "| 1 | MTAP + CDKN2A 9p21 co-deletion | PRMT5 inhibitors (MRTX1719 class) | TOY_MTAP_NULL | PRMT5 inhibitors should appear in top 10; CDK4/6 inhibitors via CDKN2A loss are also valid for this tumor and may rank higher |"
    echo "| 2 | CDK4 amplification | CDK4/6 inhibitors (palbociclib, ribociclib, abemaciclib) | TOY_FP_CDK4 | CDK4/6 inhibitors should be top 3 |"
    echo "| 3 | FGFR4 V550L hotspot | FGFR inhibitors (erdafitinib, futibatinib, FGF401-class) | TOY_TUMOR | erdafitinib + futibatinib in top 5 |"
    echo "| 4 | NRAS Q61K (with PTEN context) | MEK inhibitors (trametinib, selumetinib) | TOY_TUMOR | trametinib + selumetinib at #1-2 |"
    echo
    echo "Per-sample reports: see \`results/<sample_id>/phase5/<sample_id>.report.md\`."
} >> "$SUMMARY"

echo
echo "==========================================="
echo "Wrote multisample summary to $SUMMARY"
cat "$SUMMARY"
