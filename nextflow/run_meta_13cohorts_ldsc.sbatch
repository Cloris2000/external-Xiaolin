#!/bin/bash
#SBATCH --job-name=meta_13cohorts_ldsc
#SBATCH --partition=medium
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=logs/meta_13cohorts_ldsc_%j.out
#SBATCH --error=logs/meta_13cohorts_ldsc_%j.err

NF_DIR=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow
cd ${NF_DIR}/run_dirs/meta_13cohorts

echo "Starting 13-cohort LDSC SNP-heritability estimation..."
echo "Reuses harmonized sumstats from: ${NF_DIR}/results/meta_analysis_13cohorts/"

nextflow run ${NF_DIR}/standalone_meta_analysis.nf \
    -c ${NF_DIR}/nextflow.config.standalone_meta_13cohorts_ldsc \
    -work-dir ${NF_DIR}/work/meta_13cohorts

echo ""
echo "LDSC complete. Results: ${NF_DIR}/results/meta_analysis_13cohorts_design_matrix/ldsc/"
