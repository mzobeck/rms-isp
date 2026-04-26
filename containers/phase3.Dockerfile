# Dockerfile for phase 3 (DepMap dependency + OpenPedCan expression integration).
#
# v0.10 status: phase 3 uses pure-stdlib Python and runs from rms-isp/base:0.10.0.
# The DepMap and OpenPedCan summaries are pre-baked TSVs in assets/, so phase 3
# itself is just a TSV-join. The fetcher scripts (bin/fetch_depmap.py,
# bin/fetch_openpedcan_expression.py) need pyreadr + curl but those are
# refresh-only tools, not pipeline runtime.
#
# Anticipated v0.11+ additions if we move per-tumor expression analysis here:
#   - bcbio / Salmon for quantification from RNA-seq BAMs
#   - GSEApy or similar for pathway-level scoring

FROM rms-isp/base:0.10.0
LABEL pipeline="rms-isp"
LABEL phase="phase3_dependency"
LABEL version="0.10.0"
