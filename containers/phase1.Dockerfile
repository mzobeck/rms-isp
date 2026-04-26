# Dockerfile for phase 1 (variant + CNA + fusion annotation).
#
# v0.10 status: phase 1 uses pure-stdlib Python and runs from rms-isp/base:0.10.0.
# This per-phase Dockerfile is a placeholder for v0.11+ when phase 1 grows real
# annotation tooling. Anticipated additions when phase 1 implements the pilot
# plan §6.1 spec:
#   - Ensembl VEP 111 + cache (~50 GB; should mount as volume rather than bake in)
#   - bcftools (variant normalization)
#   - GATK4 / Mutect2 (re-calling from BAMs)
#   - cnvkit + FACETS (CNA calling from BAMs)
#   - STAR-Fusion / Arriba (fusion calling from RNA-seq BAMs)

FROM rms-isp/base:0.10.0
LABEL pipeline="rms-isp"
LABEL phase="phase1_variants"
LABEL version="0.10.0"
# v0.11+: install VEP + bcftools here.
