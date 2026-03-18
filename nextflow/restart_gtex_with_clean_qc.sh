#!/bin/bash
#
# Restart GTEx pipeline with clean genotype QC outputs
# This ensures CREATE_PGEN re-runs with sample filtering (--samples_to_keep)
#

set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

echo "=============================================="
echo "Restarting GTEx pipeline with clean QC"
echo "=============================================="
echo ""

# 1. Kill current run if it's running
echo "Step 1: Stopping current run..."
if pgrep -f "run_8_cohorts_truly_parallel" >/dev/null 2>&1; then
    pkill -f "run_8_cohorts_truly_parallel" && echo "  Killed run_8_cohorts script"
    sleep 2
fi
if pgrep -f "nextflow.*combined_pipeline_v2" >/dev/null 2>&1; then
    pkill -f "nextflow.*combined_pipeline_v2" && echo "  Killed Nextflow processes"
    sleep 2
fi
echo "  Done"
echo ""

# 2. Delete Nextflow work and cache (to force PHENO_PREP regeneration with new wgs_psam_file)
echo "Step 2: Cleaning Nextflow work and cache directories..."
ISOLATED_BASE="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/isolated_runs/gtex"
if [ -d "$ISOLATED_BASE/work" ]; then
    echo "  Removing work directory..."
    rm -rf "$ISOLATED_BASE/work"
    echo "    Removed"
fi
if [ -d "$ISOLATED_BASE/cache" ]; then
    echo "  Removing cache directory..."
    rm -rf "$ISOLATED_BASE/cache"
    echo "    Removed"
fi
echo "  Done (pipeline will run from scratch with -resume, but no cache reuse)"
echo ""

# 3. Delete old genotype QC outputs (so they get regenerated with sample filtering)
echo "Step 3: Removing old genotype QC outputs..."
GTX_QC_WORK="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/isolated_runs/gtex/launch/work/GTEx"
GTX_QC_RESULTS="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/isolated_runs/gtex/launch/results/GTEx"

if [ -d "$GTX_QC_WORK" ]; then
    echo "  Removing per-chromosome pgen files in work dir..."
    rm -f "$GTX_QC_WORK"/GTEx.QC.*.{pgen,psam,pvar,pvar.zst,log}
    rm -f "$GTX_QC_WORK"/*.prune.* "$GTX_QC_WORK"/*.het "$GTX_QC_WORK"/*.valid
    echo "    Removed files from work dir"
fi

if [ -d "$GTX_QC_RESULTS" ]; then
    echo "  Removing final QC outputs in results dir..."
    rm -f "$GTX_QC_RESULTS"/GTEx.QC.*.{pgen,psam,pvar,log}
    echo "    Removed merged/final QC files"
fi
echo "  Done"
echo ""

# 4. Restart pipeline
echo "Step 4: Restarting GTEx pipeline..."
echo "  Using: run_8_cohorts_truly_parallel.sh"
echo "  Log: logs/truly_parallel_$(date +%Y%m%d_%H%M%S)/gtex.log"
echo ""
bash run_8_cohorts_truly_parallel.sh
