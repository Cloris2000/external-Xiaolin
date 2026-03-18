#!/bin/bash
#SBATCH --job-name=meta_sens_loco_NIMH_HBCC_h650
#SBATCH --partition=medium
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/logs/meta_sensitivity/loco_NIMH_HBCC_h650_%j.out
#SBATCH --error=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/logs/meta_sensitivity/loco_NIMH_HBCC_h650_%j.err

cd /external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow
nextflow run /external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/standalone_meta_analysis.nf \
    -c /external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/sensitivity_configs/nextflow.config.sensitivity.loco_NIMH_HBCC_h650 \
    -work-dir /external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/work/meta_sensitivity/loco_NIMH_HBCC_h650 \
    -resume
