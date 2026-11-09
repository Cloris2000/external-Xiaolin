#!/bin/bash
# 5-fold out-of-fold scores for the 3 consensus models (see cons_oof.py).
# Shreejoy's checkout/env read-only; writes to /scratch/zhoux156 only.
#   sbatch --partition=debug --time=0:45:00 cons_oof.sbatch
#SBATCH --job-name=celmod_cons_oof
#SBATCH --account=def-shreejoy
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=/scratch/zhoux156/logs/celmod_strategies/%x_%j.out
set -euo pipefail
HERE=/project/rrg-shreejoy/zhoux156/external-Xiaolin/ctp_reliability_paper/celmod_strategies
export CTCOMP_REPO=/scratch/shreejoy/composition/wt_gtex
export REGISTRY_RELEASE=2026-10-21c
source "$CTCOMP_REPO/env.sh"
unset LABCORE_READS
export THREADS=16
export OMP_NUM_THREADS="$THREADS" MKL_NUM_THREADS="$THREADS"
export OPENBLAS_NUM_THREADS="$THREADS" NUMEXPR_NUM_THREADS="$THREADS"
export LEVELS=supertype BULK_SRC=light BULK_TGT=raw
export STRATEGY_OUT=/scratch/zhoux156/ctp_reliability_paper/celmod_strategies/matrix
"$PY" "$HERE/cons_oof.py"
