#!/bin/bash
#SBATCH --job-name=rosmap_array_vcfs
#SBATCH --output=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/logs/rosmap_array_vcfs_%j.out
#SBATCH --time=4:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1

set -euo pipefail

mkdir -p /external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/logs

source ~/\.anaconda3/etc/profile.d/conda.sh
conda activate bcftools_env 2>/dev/null || true

bash /external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/scripts/prep_rosmap_array_vcfs.sh
