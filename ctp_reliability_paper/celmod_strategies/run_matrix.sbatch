#!/bin/bash
# CelMod trainer x target matrix for the training-strategy figure (see run_matrix.py).
# Uses Shreejoy's composition checkout/env read-only; writes to /scratch/zhoux156 only.
#
# Usage:
#   smoke (reproduces Shreejoy's 2026-10-21c rows; his manifest, MIN_FIT 90, 5 CORE types):
#     sbatch --partition=debug --time=0:45:00 run_matrix.sbatch smoke
#   full (our manifest, MIN_FIT 89, all supertypes):
#     sbatch --time=4:00:00 run_matrix.sbatch full
#SBATCH --job-name=celmod_strategies
#SBATCH --account=def-shreejoy
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=/scratch/zhoux156/logs/celmod_strategies/%x_%j.out
set -euo pipefail
MODE="${1:?usage: run_matrix.sbatch smoke|full}"
HERE=/project/rrg-shreejoy/zhoux156/external-Xiaolin/ctp_reliability_paper/celmod_strategies
BASE=/scratch/zhoux156/ctp_reliability_paper/celmod_strategies

# the checkout and release that built results_2026-10-21c/retention (its DONE.json)
export CTCOMP_REPO=/scratch/shreejoy/composition/wt_gtex
export REGISTRY_RELEASE=2026-10-21c
source "$CTCOMP_REPO/env.sh"     # sets CTCOMP_RESULTS (read only here), CTCOMP_BULK, PY
unset LABCORE_READS

export THREADS=16                # as Shreejoy's retention_matrix.sbatch
export OMP_NUM_THREADS="$THREADS" MKL_NUM_THREADS="$THREADS"
export OPENBLAS_NUM_THREADS="$THREADS" NUMEXPR_NUM_THREADS="$THREADS"
export LEVELS=supertype BULK_SRC=light BULK_TGT=raw

case "$MODE" in
  smoke)
    export STRATEGY_MANIFEST="$CTCOMP_REPO/manifests/retention_pairs.tsv"
    export STRATEGY_OUT="$BASE/smoke" MIN_FIT=90 MAX_TYPES=5 ;;
  full)
    export STRATEGY_MANIFEST="$HERE/strategy_pairs.tsv"
    export STRATEGY_OUT="$BASE/matrix" MIN_FIT=89 MAX_TYPES=0 ;;
  *) echo "unknown mode $MODE" >&2; exit 2 ;;
esac
mkdir -p "$STRATEGY_OUT"
"$PY" "$HERE/run_matrix.py"
