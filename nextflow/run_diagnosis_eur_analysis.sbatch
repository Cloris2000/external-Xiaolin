#!/bin/bash
# Submit ancestry-matched diagnosis stratified metas + comparison.
#
# AD_EUR:           ROSMAP, ROSMAP_array, Mayo, MSBB (4 EUR AD cohorts)
# psychiatric_EUR:  CMC_MSSM, CMC_PENN, CMC_PITT (3 EUR psychiatric cohorts)
#
# Both strata run in parallel as 19-cell-type SLURM arrays; comparison runs
# after both merge jobs complete.
#
# Usage:
#   bash run_diagnosis_eur_analysis.sbatch

set -euo pipefail

NF_DIR=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow
mkdir -p "${NF_DIR}/logs/meta_sensitivity"
cd "${NF_DIR}"

echo "Submitting AD_EUR stratified meta array (19 cell types)..."
AD_ARRAY=$(sbatch --parsable "${NF_DIR}/run_stratified_meta_ad_eur_array.sbatch")
echo "  Array job ${AD_ARRAY}"

echo "Submitting psychiatric_EUR stratified meta array (19 cell types, parallel with AD_EUR)..."
PSY_ARRAY=$(sbatch --parsable "${NF_DIR}/run_stratified_meta_psychiatric_eur_array.sbatch")
echo "  Array job ${PSY_ARRAY}"

echo "Submitting AD_EUR merge (after array ${AD_ARRAY})..."
AD_MERGE=$(sbatch --parsable --dependency=afterok:${AD_ARRAY} "${NF_DIR}/run_stratified_meta_ad_eur_merge.sbatch")
echo "  Merge job ${AD_MERGE}"

echo "Submitting psychiatric_EUR merge (after array ${PSY_ARRAY})..."
PSY_MERGE=$(sbatch --parsable --dependency=afterok:${PSY_ARRAY} "${NF_DIR}/run_stratified_meta_psychiatric_eur_merge.sbatch")
echo "  Merge job ${PSY_MERGE}"

echo "Submitting ancestry-matched diagnosis comparison (after both merges)..."
CMP_JOB=$(sbatch --parsable --dependency=afterok:${AD_MERGE}:${PSY_MERGE} "${NF_DIR}/run_compare_stratified_meta_diagnosis_eur.sbatch")
echo "  Comparison job ${CMP_JOB}"

echo ""
echo "Monitor: squeue -u \$USER | grep -E 'meta_AD_EUR|meta_psych|cmp_diag'"
echo "Results: results/meta_sensitivity/stratified_comparison_diagnosis_eur/"
