#!/bin/bash
#SBATCH --job-name=sn_mssm
#SBATCH --partition=medium
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=logs/sn_mssm_%j.out
#SBATCH --error=logs/sn_mssm_%j.err

NF_DIR=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow
cd ${NF_DIR}/run_dirs/sn_CMC_MSSM

export JAVA_HOME="$HOME/.sdkman/candidates/java/current"
export PATH="$HOME/.local/bin:$JAVA_HOME/bin:$PATH"

nextflow run ${NF_DIR}/combined_pipeline_v2.nf \
  -c ${NF_DIR}/nextflow.config.combined.sn_mssm \
  -w ${NF_DIR}/work/sn_CMC_MSSM \
  -resume
