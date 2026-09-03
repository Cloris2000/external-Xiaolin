#!/bin/bash
#SBATCH --job-name=sn_hbcc_Omni5M
#SBATCH --partition=compute
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=logs/sn_hbcc_Omni5M_%j.out
#SBATCH --error=logs/sn_hbcc_Omni5M_%j.err
#SBATCH --account=rrg-shreejoy

source "${SITE_ENV:-/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/site_env.sh}"
cd ${NF_DIR}/run_dirs/sn_NIMH_HBCC_Omni5M

export JAVA_HOME="$HOME/.sdkman/candidates/java/current"
export PATH="$HOME/.local/bin:$JAVA_HOME/bin:$PATH"

nextflow run ${NF_DIR}/combined_pipeline_v2.nf \
  -c ${NF_DIR}/nextflow.config.combined.sn_hbcc_Omni5M \
  -w ${NF_DIR}/work/sn_NIMH_HBCC_Omni5M \
  -resume
