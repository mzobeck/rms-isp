# Dockerfile for phase 5 (confidence scoring + Markdown report rendering).
#
# v0.10 status: phase 5 uses pure-stdlib Python and runs from rms-isp/base:0.10.0.
# Anticipated v0.11+ additions:
#   - WeasyPrint or similar for PDF report rendering
#   - Jinja2 templating if Markdown layout outgrows the current f-string approach
#   - PDF signing for clinical traceability (v1.0+)

FROM rms-isp/base:0.10.0
LABEL pipeline="rms-isp"
LABEL phase="phase5_scoring"
LABEL version="0.10.0"
