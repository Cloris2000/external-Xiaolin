#!/usr/bin/env bash
# Copy pipeline results from scratch to the durable /project tree.
#
# Trillium compute nodes mount /project read-only, so jobs publish to
# $RESULTS_ROOT on scratch.  /scratch is purged periodically, so results are not
# durable until they have been copied across.  This has to run on a LOGIN node,
# which is the only place /project is writable.
#
#   scripts/sync_results_to_project.sh                       # everything
#   scripts/sync_results_to_project.sh meta_analysis_15cohorts
#   DRY_RUN=1 scripts/sync_results_to_project.sh meta_analysis_15cohorts

set -euo pipefail

source "${SITE_ENV:-/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/site_env.sh}"

if [ ! -w "${NF_DIR}" ]; then
    echo "ERROR: ${NF_DIR} is not writable from $(hostname)." >&2
    echo "This must run on a login node; /project is read-only on compute nodes." >&2
    exit 1
fi

subdirs=("$@")
if [ ${#subdirs[@]} -eq 0 ]; then
    mapfile -t subdirs < <(cd "${RESULTS_ROOT}" 2>/dev/null && ls -1)
fi

if [ ${#subdirs[@]} -eq 0 ]; then
    echo "Nothing to sync: ${RESULTS_ROOT} is empty."
    exit 0
fi

rsync_opts=(-a --human-readable --info=stats2)
[ -n "${DRY_RUN:-}" ] && rsync_opts+=(--dry-run) && echo "(dry run)"

for sub in "${subdirs[@]}"; do
    src="${RESULTS_ROOT}/${sub}"
    dst="${PROJECT_RESULTS}/${sub}"
    if [ ! -d "${src}" ]; then
        echo "SKIP  ${sub}: not present under ${RESULTS_ROOT}" >&2
        continue
    fi
    echo "==> ${src}  ->  ${dst}"
    mkdir -p "${dst}"
    # Trailing slash on src copies the contents, not the directory itself.
    rsync "${rsync_opts[@]}" "${src}/" "${dst}/"
done

echo "Sync complete."
