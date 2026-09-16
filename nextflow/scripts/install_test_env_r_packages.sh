#!/usr/bin/env bash
# Add the three R packages the imported `test` env does not ship.
#
# Use mamba, not conda. conda's classic solver can sit for hours against this
# env (and is what hung the original install). mamba resolves the same
# conda-forge packages without pulling the rest of the env apart.
#
#   r-coloc and r-susier  — conda-forge
#   topr                  — not packaged as r-topr on conda-forge/bioconda,
#                           so it is installed from CRAN into the same env
#
# Safe to re-run: already-present packages are left alone.
# Must be run on a login node (compute nodes have no outbound network).

set -euo pipefail

CONDA_ENV="${CONDA_ENV:-test}"
CONDA_ROOT="${CONDA_ROOT:-$HOME/miniforge3}"

# shellcheck disable=SC1091
source "${CONDA_ROOT}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV}"

if ! command -v mamba >/dev/null 2>&1; then
    echo "ERROR: mamba is not on PATH. Install mamba in the base env first." >&2
    exit 1
fi

need_conda=()
Rscript -e 'quit(status = as.integer(!requireNamespace("coloc",  quietly=TRUE)))' || need_conda+=(r-coloc)
Rscript -e 'quit(status = as.integer(!requireNamespace("susieR", quietly=TRUE)))' || need_conda+=(r-susier)

if [ "${#need_conda[@]}" -gt 0 ]; then
    echo "mamba install -n ${CONDA_ENV} -c conda-forge ${need_conda[*]}"
    mamba install -n "${CONDA_ENV}" -y -c conda-forge "${need_conda[@]}"
else
    echo "r-coloc and r-susier already present; skipping mamba"
fi

if ! Rscript -e 'quit(status = as.integer(!requireNamespace("topr", quietly=TRUE)))'; then
    echo "Installing topr from CRAN (no r-topr package on conda-forge)"
    Rscript -e 'install.packages("topr", repos="https://cloud.r-project.org")'
else
    echo "topr already present; skipping CRAN"
fi

Rscript -e '
cat(sprintf("coloc  %s\n", as.character(packageVersion("coloc"))))
cat(sprintf("susieR %s\n", as.character(packageVersion("susieR"))))
cat(sprintf("topr   %s\n", as.character(packageVersion("topr"))))
'
