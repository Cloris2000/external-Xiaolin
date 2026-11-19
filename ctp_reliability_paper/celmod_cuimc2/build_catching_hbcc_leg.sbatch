#!/bin/bash
# CelMod out-of-fold leg NIMH_HBCC bulk x Catching HBCC, EXPLORATORY MIN_PAIRED=80 (see build_catching_hbcc_leg.py).
# Uses Shreejoy's composition checkout and environment read-only; writes to scratch only.
#SBATCH --job-name=catching_hbcc_celmod
#SBATCH --account=def-shreejoy
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=0:45:00
#SBATCH --output=/scratch/zhoux156/logs/celmod_catching_hbcc/%x_%j.out
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

export CATCHING_HBCC_LEG_OUT=/scratch/zhoux156/ctp_reliability_paper/celmod_catching_hbcc
mkdir -p "$CATCHING_HBCC_LEG_OUT"
"$PY" /project/rrg-shreejoy/zhoux156/external-Xiaolin/ctp_reliability_paper/celmod_cuimc2/build_catching_hbcc_leg.py
