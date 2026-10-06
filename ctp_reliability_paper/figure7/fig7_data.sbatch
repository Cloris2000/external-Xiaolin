#!/bin/bash
#SBATCH --job-name=fig7_data
#SBATCH --account=def-shreejoy
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --output=/scratch/zhoux156/logs/fig7_all_supertypes/fig7_data_%j.out
set -euo pipefail
# /project is read-only on Trillium compute nodes (nextflow/site_env.sh), so
# output must go to scratch; copy back to figure7/data/ from a login node.
export FIG7_OUT=/scratch/zhoux156/ctp_reliability_paper/figure7_data
mkdir -p "$FIG7_OUT"
cd /project/rrg-shreejoy/zhoux156/external-Xiaolin/ctp_reliability_paper/figure7
PYTHONNOUSERSITE=1 /home/zhoux156/miniforge3/envs/brisc/bin/python fig7_data.py
