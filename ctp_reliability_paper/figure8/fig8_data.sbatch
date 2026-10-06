#!/bin/bash
#SBATCH --job-name=fig8_data
#SBATCH --account=def-shreejoy
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --output=/scratch/zhoux156/logs/fig8_depth/fig8_data_%j.out
set -euo pipefail
# /project is read-only on Trillium compute nodes (nextflow/site_env.sh), so
# output must go to scratch; copy back to figure8/data/ from a login node.
export FIG8_OUT=/scratch/zhoux156/ctp_reliability_paper/figure8_data
mkdir -p "$FIG8_OUT"
cd /project/rrg-shreejoy/zhoux156/external-Xiaolin/ctp_reliability_paper/figure8
PYTHONNOUSERSITE=1 /home/zhoux156/miniforge3/envs/brisc/bin/python fig8_data.py
