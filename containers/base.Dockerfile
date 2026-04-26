# RMS-ISP base image: python:3.11-slim plus procps (Nextflow needs `ps` for
# task-metric collection). All five phase processes share this image at v0.10
# because the pipeline is pure-stdlib Python. Phase-specific overrides start
# in v0.11+ when phase 1 gets VEP, phase 2 gets Boltz/Chai, etc.
#
# Build:
#   docker build -t rms-isp/base:0.10.0 -f containers/base.Dockerfile .
# Reference from Nextflow:
#   process.container = 'rms-isp/base:0.10.0'

FROM python:3.11-slim

LABEL maintainer="Mark Zobeck <mzobeck@gmail.com>"
LABEL pipeline="rms-isp"
LABEL stage="base"
LABEL version="0.10.0"

RUN apt-get update && apt-get install -y --no-install-recommends \
    procps ca-certificates && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /work
CMD ["python3"]
