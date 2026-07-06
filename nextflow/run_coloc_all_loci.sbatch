#!/bin/bash
#SBATCH --job-name=coloc_all_loci
#SBATCH --array=1-19%5
#SBATCH --time=6:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/coloc/coloc_full_%A_%a.out
#SBATCH --error=logs/coloc/coloc_full_%A_%a.err

# Colocalize ALL suggestive loci (p < 1e-5) using full 500 kb windows.
# Differences from run_coloc.sbatch:
#   - All SNPs in each locus window (not just suggestive SNPs)
#   - min_snps = 10 (coloc.abf practical minimum; 50 was overly conservative)
#
# Usage:
#   sbatch run_coloc_all_loci.sbatch
#
# After all jobs finish:
#   Rscript scripts/coloc/04_summarize_coloc.R \
#     --results_dir results/coloc/coloc_results_full \
#     --loci_dir    results/coloc/loci_full \
#     --output_dir  results/coloc/full

set -euo pipefail

export R_LIBS="/external/rprshnas01/kcni/xzhou/.anaconda3/lib/R/library"
export PATH="/nethome/kcni/xzhou/.anaconda3/bin:$PATH"

PROJECT_DIR="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow"
RESULTS_DIR="$PROJECT_DIR/results/meta_analysis_15cohorts"
COLOC_LOCI_DIR="$PROJECT_DIR/results/coloc/loci_full"
COLOC_RES_DIR="$PROJECT_DIR/results/coloc/coloc_results_full"
DISEASE_DIR="$PROJECT_DIR/results/coloc/disease_gwas"
SCRIPTS_DIR="$PROJECT_DIR/scripts/coloc"
CELL_TYPE_N=2812
MIN_SNPS=10

mkdir -p logs/coloc "$COLOC_LOCI_DIR" "$COLOC_RES_DIR"

CELL_TYPES=(
  "Astrocyte" "Endothelial" "IT" "L4.IT" "L5.6.IT.Car3" "L5.6.NP"
  "L5.ET" "L6.CT" "L6b" "LAMP5" "Microglia" "OPC" "Oligodendrocyte"
  "PAX6" "PVALB" "Pericyte" "SST" "VIP" "VLMC"
)

IDX=$(( SLURM_ARRAY_TASK_ID - 1 ))
CELL_TYPE="${CELL_TYPES[$IDX]}"

echo "============================================"
echo "Job:       $SLURM_JOB_ID (array task $SLURM_ARRAY_TASK_ID)"
echo "Cell type: $CELL_TYPE"
echo "Node:      $(hostname)"
echo "Started:   $(date)"
echo "============================================"

TBL_FILE=$(ls "$RESULTS_DIR"/"${CELL_TYPE}"_meta_analysis_*.tbl 2>/dev/null | \
           grep -v "[0-9]\.tbl$" | head -1)
if [ -z "$TBL_FILE" ]; then
  echo "ERROR: No .tbl file found for $CELL_TYPE"
  exit 1
fi
echo "TBL file: $TBL_FILE"

LOCI_FILE="$COLOC_LOCI_DIR/${CELL_TYPE}_loci.tsv"
LOCUS_DATA="$COLOC_LOCI_DIR/${CELL_TYPE}_locus_data.tsv.gz"
COLOC_OUT="$COLOC_RES_DIR/${CELL_TYPE}_coloc_results.tsv"

echo "--- Step 1: Extract loci (full-window SNPs for coloc) ---"
Rscript "$SCRIPTS_DIR/01_extract_loci.R" \
  --cell_type         "$CELL_TYPE" \
  --tbl_file          "$TBL_FILE" \
  --output_dir        "$COLOC_LOCI_DIR" \
  --p_threshold       1e-5 \
  --window_kb         500 \
  --full_window_snps  TRUE

NLOCI=$(tail -n +2 "$LOCI_FILE" | wc -l)
echo "  $NLOCI loci; locus data: $LOCUS_DATA"

echo "--- Step 2: coloc.abf (min_snps=$MIN_SNPS) ---"
Rscript "$SCRIPTS_DIR/03_run_coloc.R" \
  --cell_type    "$CELL_TYPE" \
  --loci_file    "$LOCI_FILE" \
  --locus_data   "$LOCUS_DATA" \
  --disease_dir  "$DISEASE_DIR" \
  --output_dir   "$COLOC_RES_DIR" \
  --cell_type_N  "$CELL_TYPE_N" \
  --run_susie    FALSE \
  --min_snps     "$MIN_SNPS"

echo "============================================"
echo "Finished: $(date)"
echo "Output:   $COLOC_OUT"
echo "============================================"
