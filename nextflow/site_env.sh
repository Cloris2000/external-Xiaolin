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

# Nextflow work dirs and other scratch output.  /scratch is purged periodically,
# so nothing here may be treated as a durable result.
WORK_ROOT="${WORK_ROOT:-/scratch/${USER}/nf_work}"

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

export NF_DIR SCC_DIR SCC_ROOT REFS_ROOT DATA_ROOT WORK_ROOT
export SBATCH_ACCOUNT SBATCH_PARTITION NODE_CPUS
export CONDA_ENV CONDA_ROOT RSCRIPT PYTHON NEXTFLOW_MODULE
