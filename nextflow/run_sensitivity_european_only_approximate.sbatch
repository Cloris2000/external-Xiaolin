#!/bin/bash
#SBATCH --job-name=meta_sens_european_only_approximate
#SBATCH --partition=medium
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/logs/meta_sensitivity/european_only_approximate_%j.out
#SBATCH --error=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/logs/meta_sensitivity/european_only_approximate_%j.err

cd /external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow
nextflow run /external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/standalone_meta_analysis.nf \
    -c /external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/sensitivity_configs/nextflow.config.sensitivity.european_only_approximate \
    -work-dir /external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/work/meta_sensitivity/european_only_approximate \
    -resume
