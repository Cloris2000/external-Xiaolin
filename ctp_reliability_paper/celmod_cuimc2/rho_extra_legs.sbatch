#!/bin/bash
# Spearman rho for the extra Fig 4 legs from saved predictions (see rho_extra_legs.py); no CelMod fit.
# Uses Shreejoy's composition checkout and environment read-only; writes to scratch only.
#SBATCH --job-name=rho_extra_legs
#SBATCH --account=def-shreejoy
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=0:20:00
#SBATCH --output=/scratch/zhoux156/logs/celmod_rho/%x_%j.out
set -euo pipefail
# the checkout and release that built results_2026-10-21c/deconv (DONE.json)
export CTCOMP_REPO=/scratch/shreejoy/composition/wt_gtex
export REGISTRY_RELEASE=2026-10-21c
source "$CTCOMP_REPO/env.sh"     # sets CTCOMP_RESULTS (read), CTCOMP_BULK, PY, PYTHONPATH
unset LABCORE_READS              # no read-tracking files in Shreejoy's results root

# same thread setting as build_legs.sbatch
export THREADS=16
export OMP_NUM_THREADS="$THREADS" MKL_NUM_THREADS="$THREADS"
export OPENBLAS_NUM_THREADS="$THREADS" NUMEXPR_NUM_THREADS="$THREADS"

export CUIMC2_LEG_OUT=/scratch/zhoux156/ctp_reliability_paper/celmod_cuimc2
mkdir -p "$CUIMC2_LEG_OUT"
"$PY" /project/rrg-shreejoy/zhoux156/external-Xiaolin/ctp_reliability_paper/celmod_cuimc2/rho_extra_legs.py
