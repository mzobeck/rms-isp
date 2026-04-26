# Dockerfile for phase 4 (drug matching).
#
# v0.10 status: phase 4 uses pure-stdlib Python and runs from rms-isp/base:0.10.0.
# The DGIdb and ClinicalTrials.gov caches are pre-baked TSVs; phase 4 itself is
# a TSV-join over (gene, drug). Fetchers (bin/fetch_dgidb.py,
# bin/fetch_clinicaltrials.py) only need urllib from stdlib so they work
# in this image too.
#
# Anticipated v0.11+ additions:
#   - signature-based matching via CMap / LINCS L1000 client
#   - OpenTargets GraphQL fetcher (currently skipped due to v4 schema complexity)
#   - Molecular Targets Platform integration

FROM rms-isp/base:0.10.0
LABEL pipeline="rms-isp"
LABEL phase="phase4_drugs"
LABEL version="0.10.0"
