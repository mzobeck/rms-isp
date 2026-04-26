#!/usr/bin/env bash
# Build the RMS-ISP container images locally.
#
# v0.10: only the base image is functionally distinct (stdlib Python + procps).
# The five per-phase Dockerfiles all extend the base; building them in v0.10
# produces images that are byte-identical to the base, but the per-phase tags
# document the v0.11+ overlay points (VEP for phase 1, Boltz/Chai for phase 2,
# etc.).
#
# Usage:
#   containers/build.sh            # build base only (sufficient for v0.10)
#   containers/build.sh --all      # also build the 5 per-phase tags

set -euo pipefail

cd "$(dirname "$0")/.."

VERSION=$(awk -F"'" '/version/ { print $2; exit }' nextflow.config)
echo "==> building rms-isp/base:$VERSION"
docker build -t "rms-isp/base:$VERSION" -f containers/base.Dockerfile .

if [[ "${1:-}" == "--all" ]]; then
    for phase in phase1 phase2 phase3 phase4 phase5; do
        echo "==> building rms-isp/$phase:$VERSION"
        docker build -t "rms-isp/$phase:$VERSION" -f "containers/$phase.Dockerfile" .
    done
fi

echo
echo "Built images:"
docker images "rms-isp/*" --format "table {{.Repository}}:{{.Tag}}\t{{.Size}}\t{{.CreatedSince}}"
