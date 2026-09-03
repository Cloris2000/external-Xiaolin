#!/bin/bash
# Trillium site environment.
#
# Single source of truth for cluster paths, the conda environment and SLURM
# settings used by the sbatch drivers and shell wrappers in this repo.  The
# Nextflow equivalent is nextflow.config.trillium; keep the two in sync.
#
# Usage, near the top of a driver:
#     source "$(dirname "${BASH_SOURCE[0]}")/site_env.sh"
#
# Every variable honours a pre-existing value, so a caller can override any of
# them in the environment without editing this file:
#     NF_DIR=/some/other/checkout sbatch run_coloc_all_loci.sbatch

# Pipeline checkout (was /external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow on the SCC).
NF_DIR="${NF_DIR:-/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow}"

# Legacy SCC tree copied to Trillium: per-cohort regenie_step2 outputs, disease
# GWAS, sn meta tables and the LDSC European reference bundle all still live here.
SCC_DIR="${SCC_DIR:-/project/rrg-shreejoy/zhoux156/Xiaolin/SCC/nextflow}"
SCC_ROOT="${SCC_ROOT:-/project/rrg-shreejoy/zhoux156/Xiaolin/SCC}"

# Shared tools and reference data.
REFS_ROOT="${REFS_ROOT:-/project/rrg-shreejoy/pipeline_refs}"
DATA_ROOT="${DATA_ROOT:-/project/rrg-shreejoy}"

# Trillium compute nodes mount /project AND $HOME read-only; only /scratch is
# writable (verified on tri0849).  So every path a running job writes to - job
# logs, Nextflow's work and home dirs, R/tmp scratch and all pipeline results -
# has to live under SCRATCH_ROOT.  Durable results are synced back to /project
# from a login node afterwards; see scripts/sync_results_to_project.sh.
#
# /scratch is also purged periodically, so nothing under it is a durable result
# until it has been synced.
SCRATCH_ROOT="${SCRATCH_ROOT:-/scratch/${USER}}"
WORK_ROOT="${WORK_ROOT:-${SCRATCH_ROOT}/nf_work}"
RESULTS_ROOT="${RESULTS_ROOT:-${SCRATCH_ROOT}/results}"
LOG_ROOT="${LOG_ROOT:-${SCRATCH_ROOT}/logs}"

# Durable copy of the results, on read-only-at-runtime /project.
PROJECT_RESULTS="${PROJECT_RESULTS:-${NF_DIR}/results}"

# Nextflow defaults NXF_HOME to $HOME/.nextflow and writes plugins and caches
# there, which fails on a compute node.
export NXF_HOME="${NXF_HOME:-${SCRATCH_ROOT}/.nextflow}"
export TMPDIR="${TMPDIR:-${SCRATCH_ROOT}/tmp}"
mkdir -p "${WORK_ROOT}" "${RESULTS_ROOT}" "${LOG_ROOT}" "${NXF_HOME}" "${TMPDIR}" 2>/dev/null || true

# SLURM.  Trillium is SelectType=select/linear, so every job is allocated a whole
# 192-core / 767 GB node on the `compute` partition (24 h max).  The SCC partition
# names (short/medium/mediumtmp/long) do not exist here.
SBATCH_ACCOUNT="${SBATCH_ACCOUNT:-rrg-shreejoy}"
SBATCH_PARTITION="${SBATCH_PARTITION:-compute}"
NODE_CPUS="${NODE_CPUS:-192}"

# Conda environment imported from the SCC.  Supplies R 4.4.3 (data.table, dplyr,
# tidyr, optparse, ggplot2, cowplot, qqman, coloc, susieR, topr) and Python 3.11.
# Do not switch this to r_env, which is a separate, less complete environment.
CONDA_ENV="${CONDA_ENV:-test}"
CONDA_ROOT="${CONDA_ROOT:-$HOME/miniforge3}"

# `conda activate` trips over `set -u`, which most drivers enable.
if [ "$(basename "${CONDA_PREFIX:-none}")" != "${CONDA_ENV}" ]; then
    _site_env_had_nounset=0
    case "$-" in *u*) _site_env_had_nounset=1; set +u ;; esac
    # shellcheck disable=SC1091
    source "${CONDA_ROOT}/etc/profile.d/conda.sh"
    conda activate "${CONDA_ENV}"
    [ "${_site_env_had_nounset}" = "1" ] && set -u
    unset _site_env_had_nounset
fi

# The SCC drivers pointed R_LIBS at a standalone Anaconda tree that does not
# exist here.  Leaving it unset makes R resolve against the active conda env.
unset R_LIBS R_LIBS_USER 2>/dev/null || true

RSCRIPT="${RSCRIPT:-${CONDA_PREFIX}/bin/Rscript}"
PYTHON="${PYTHON:-${CONDA_PREFIX}/bin/python}"

# Nextflow 26's config parser rejects the cross-referencing style these configs
# use (e.g. output_dir = "${base_dir}/..."), so the version is pinned.
NEXTFLOW_MODULE="${NEXTFLOW_MODULE:-nextflow/25.10.2}"

# Compute nodes have no outbound network, so Nextflow's "is there a newer
# version?" check burns ~130 s on a curl timeout before every run.
export NXF_DISABLE_CHECK_LATEST=true
export CAPSULE_LOG=none

export NF_DIR SCC_DIR SCC_ROOT REFS_ROOT DATA_ROOT
export SCRATCH_ROOT WORK_ROOT RESULTS_ROOT LOG_ROOT PROJECT_RESULTS
export SBATCH_ACCOUNT SBATCH_PARTITION NODE_CPUS
export CONDA_ENV CONDA_ROOT RSCRIPT PYTHON NEXTFLOW_MODULE
