#!/bin/bash
# Submit Tier 1 heterogeneity-focused ancestry screen (19 parallel cell types).
#
# Prerequisites:
#   - stratified_hit_overlap.tsv exists (from compare_stratified_meta)
#   - EUR + HBCC stratified metas complete
#
# Usage:
#   bash run_ancestry_tier1_screening.sbatch
#
# After array completes:
#   sbatch run_ancestry_tier1_merge.sbatch
#   # or: sbatch --dependency=afterok:<ARRAY_JOB_ID> run_ancestry_tier1_merge.sbatch

set -euo pipefail

NF_DIR=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow
mkdir -p "${NF_DIR}/logs/meta_sensitivity"
cd "${NF_DIR}"

OVERLAP="results/meta_sensitivity/stratified_comparison_tier1/stratified_hit_overlap.tsv"
if [[ ! -f "${OVERLAP}" ]]; then
  echo "Missing ${OVERLAP}. Running compare_stratified_meta first..."
  python3 - <<'PY'
import csv
from pathlib import Path
strata_tsv = Path("results/meta_sensitivity/ancestry_tier1_strata.tsv")
strata_tsv.parent.mkdir(parents=True, exist_ok=True)
rows = [
    {"stratum_label": "pooled_15cohorts", "results_dir": "results/meta_analysis_15cohorts"},
    {"stratum_label": "EUR_homogeneous", "results_dir": "results/meta_sensitivity/stratified_eur_homogeneous_15cohorts"},
    {"stratum_label": "HBCC_AFR_enriched", "results_dir": "results/meta_sensitivity/stratified_hbcc_15cohorts"},
]
with open(strata_tsv, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["stratum_label", "results_dir"], delimiter="\t")
    w.writeheader()
    w.writerows(rows)
PY
  python3 scripts/compare_stratified_meta.py \
      --project-dir "${NF_DIR}" \
      --strata-tsv results/meta_sensitivity/ancestry_tier1_strata.tsv \
      --output-dir results/meta_sensitivity/stratified_comparison_tier1
fi

echo "Submitting 19 parallel cell-type array jobs..."
ARRAY_JOB=$(sbatch --parsable "${NF_DIR}/run_ancestry_tier1_array.sbatch")
echo "Array job submitted: ${ARRAY_JOB}"

echo "Submitting merge job (runs after array completes)..."
MERGE_JOB=$(sbatch --parsable --dependency=afterok:${ARRAY_JOB} "${NF_DIR}/run_ancestry_tier1_merge.sbatch")
echo "Merge job submitted: ${MERGE_JOB} (depends on ${ARRAY_JOB})"

echo ""
echo "Monitor: squeue -u \$USER | grep ancestry"
echo "Per-cell-type outputs: results/meta_sensitivity/ancestry_tier1/per_celltype/"
echo "Final outputs after merge: results/meta_sensitivity/ancestry_tier1/ancestry_masking_candidates.tsv"
