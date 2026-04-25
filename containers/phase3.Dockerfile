# Dockerfile for phase3 module.
# Skeleton placeholder. Pin tools and reference databases at the phase build step.

FROM ubuntu:22.04
LABEL maintainer="Mark Zobeck <mzobeck@gmail.com>"
LABEL pipeline="rms-isp"
LABEL phase="phase3"
LABEL version="0.0.1-scaffold"

RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates curl && \
    rm -rf /var/lib/apt/lists/*

CMD ["/bin/bash"]
