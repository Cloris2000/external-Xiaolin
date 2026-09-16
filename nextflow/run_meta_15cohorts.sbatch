#!/bin/bash
#SBATCH --job-name=meta_15cohorts
#SBATCH --account=rrg-shreejoy
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=192
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/zhoux156/logs/meta_15cohorts/meta_15cohorts_%j.out
#SBATCH --error=/scratch/zhoux156/logs/meta_15cohorts/meta_15cohorts_%j.err

# Trillium allocates whole nodes (SelectType=select/linear), so this job owns all
# 192 cores and 767 GB and Nextflow schedules tasks locally inside the allocation.
# Do not add --mem; node memory is allocated in full.
#
# The geometry has to be 1 task x 192 cpus, not 192 tasks: the batch step is
# cgroup-limited to its own task's cores, so --ntasks-per-node=192 leaves this
# script (and therefore Nextflow's local executor) pinned to a single core.
# Nextflow then reports "capacity=1" and runs every cell type sequentially.

set -euo pipefail

# NF_DIR, SCC_DIR, WORK_ROOT, RESULTS_ROOT, NXF_HOME, the `test` conda env
# (R + python for the plot stages) and NEXTFLOW_MODULE all come from here.
source "${SITE_ENV:-/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/site_env.sh}"
WORK_DIR="${WORK_ROOT}/meta_15cohorts"

# Compute nodes mount NF_DIR read-only, so Nextflow's run directory (.nextflow,
# .nextflow.log, reports) and all published output have to live on scratch.
RUN_DIR="${SCRATCH_ROOT}/run_dirs/meta_15cohorts"
OUT_DIR="${RESULTS_ROOT}/meta_analysis_15cohorts"

# Pinned to 25.10.2: Nextflow 26's config parser is stricter and rejects the
# cross-reference style these configs rely on (e.g. output_dir = "${base_dir}/...").
module load "${NEXTFLOW_MODULE}"

mkdir -p "${RUN_DIR}" "${WORK_DIR}" "${OUT_DIR}"
cd "${RUN_DIR}"

# Stage 0 reads {base_dir}/results/{cohort}/regenie_step2 from the new GWAS run.
# Prefer scratch (where compute-node jobs publish); fall back to the repo copy.
echo "Checking new GWAS results for all 15 cohorts..."
all_good=true
for cohort in ROSMAP ROSMAP_array Mayo MSBB CMC_MSSM CMC_PENN CMC_PITT GTEx_v10 NABEC NIMH_HBCC_1M NIMH_HBCC_h650 NIMH_HBCC_Omni5M GVEX AMP_AD_Rush AMP_AD_Mayo; do
    count=$(ls ${RESULTS_ROOT}/${cohort}/regenie_step2/*.regenie.raw_p 2>/dev/null | wc -l)
    src="scratch"
    if [ "$count" -ne 19 ]; then
        count=$(ls ${NF_DIR}/results/${cohort}/regenie_step2/*.regenie.raw_p 2>/dev/null | wc -l)
        src="project"
    fi
    if [ "$count" -eq 19 ]; then
        echo "  OK  ${cohort}: ${count}/19 (${src})"
    else
        echo "  FAIL ${cohort}: ${count}/19 (INCOMPLETE)"
        all_good=false
    fi
done

if [ "$all_good" = false ]; then
    echo ""
    echo "ERROR: Some cohorts have incomplete GWAS results. Aborting."
    echo "Run each cohort GWAS first; do not reuse the old SCC/pre-harmonized inputs."
    exit 1
fi

echo ""
echo "All inputs present. Starting 15-cohort meta-analysis..."

nextflow run ${NF_DIR}/standalone_meta_analysis.nf \
    -c ${NF_DIR}/nextflow.config.standalone_meta_15cohorts \
    -work-dir ${WORK_DIR} \
    --output_dir ${OUT_DIR} \
    --base_dir ${NF_DIR}

echo ""
echo "15-cohort meta-analysis complete."
echo "Results written to: ${OUT_DIR}/"
echo "Sync them to ${PROJECT_RESULTS}/ from a login node with:"
echo "  ${NF_DIR}/scripts/sync_results_to_project.sh meta_analysis_15cohorts"
