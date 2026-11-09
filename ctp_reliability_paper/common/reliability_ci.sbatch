#!/bin/bash
#SBATCH --job-name=reliability_ci
#SBATCH --account=def-shreejoy
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --output=/scratch/zhoux156/logs/reliability_ci/reliability_ci_%j.out
set -euo pipefail
# /project is read-only on Trillium compute nodes: write to scratch, copy back from a login node.
export RELIABILITY_CI_OUT=/scratch/zhoux156/ctp_reliability_paper/reliability_ci
mkdir -p "$RELIABILITY_CI_OUT"
cd /project/rrg-shreejoy/zhoux156/external-Xiaolin/ctp_reliability_paper/common
PYTHONNOUSERSITE=1 /home/zhoux156/miniforge3/envs/brisc/bin/python reliability_ci.py
