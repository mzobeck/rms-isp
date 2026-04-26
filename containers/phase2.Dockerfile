# Dockerfile for phase 2 (structural reference + future docking).
#
# v0.10 status: phase 2 uses pure-stdlib Python and runs from rms-isp/base:0.10.0.
# Anticipated v0.11+ additions per pilot plan §6.2:
#   - Boltz-1 (MIT) for mutant-protein structure prediction
#   - Chai-1 (Chai Discovery) as the alternate predictor
#   - AlphaFold-Multimer for fusion junction structure
#   - AutoDock Vina for ligand docking
# All of these are GPU-bound; this Dockerfile will likely become a CUDA base
# image (nvidia/cuda + python + Boltz/Chai weights).

FROM rms-isp/base:0.10.0
LABEL pipeline="rms-isp"
LABEL phase="phase2_structure"
LABEL version="0.10.0"
# v0.11+: switch base to nvidia/cuda and install Boltz / Chai / AutoDock here.
