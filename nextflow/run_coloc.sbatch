#!/bin/bash
#SBATCH --job-name=coloc_15cohorts
#SBATCH --array=1-19%5            # 19 cell types, max 5 concurrent
#SBATCH --time=4:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/coloc/coloc_%A_%a.out
#SBATCH --error=logs/coloc/coloc_%A_%a.err

# ============================================================================
# Colocalization pipeline: 15-cohort cell-type GWAS × disease GWAS
# Runs 03_run_coloc.R for each of the 19 cell types as a SLURM array job.
#
# Prerequisites:
#   1. Run 01_extract_loci.R for all cell types (or let this script do it)
#   2. Run 02_download_disease_gwas.sh to download disease GWAS sumstats
#   3. Install R packages: coloc, susieR, data.table, dplyr, tidyr, optparse
#
# Usage:
#   sbatch run_coloc.sbatch
#
# After all array jobs complete, run the summarization:
#   Rscript scripts/coloc/04_summarize_coloc.R
# ============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# R library path — point to the conda env where coloc/susieR are installed
# ---------------------------------------------------------------------------
export R_LIBS="/external/rprshnas01/kcni/xzhou/.anaconda3/lib/R/library"
export PATH="/nethome/kcni/xzhou/.anaconda3/bin:$PATH"

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT_DIR="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow"
RESULTS_DIR="$PROJECT_DIR/results/meta_analysis_15cohorts"
COLOC_LOCI_DIR="$PROJECT_DIR/results/coloc/loci"
COLOC_RES_DIR="$PROJECT_DIR/results/coloc/coloc_results"
DISEASE_DIR="$PROJECT_DIR/results/coloc/disease_gwas"
SCRIPTS_DIR="$PROJECT_DIR/scripts/coloc"

# Total effective meta-analysis N (sum of final_gwas_n across 12 cohorts with data)
CELL_TYPE_N=2812

mkdir -p logs/coloc "$COLOC_LOCI_DIR" "$COLOC_RES_DIR"

# ---------------------------------------------------------------------------
# Cell type list (19 cell types, order matches SLURM array index 1-19)
# ---------------------------------------------------------------------------
CELL_TYPES=(
  "Astrocyte"
  "Endothelial"
  "IT"
  "L4.IT"
  "L5.6.IT.Car3"
  "L5.6.NP"
  "L5.ET"
  "L6.CT"
  "L6b"
  "LAMP5"
  "Microglia"
  "OPC"
  "Oligodendrocyte"
  "PAX6"
  "PVALB"
  "Pericyte"
  "SST"
  "VIP"
  "VLMC"
)

# SLURM_ARRAY_TASK_ID is 1-based; array index is 0-based
IDX=$(( SLURM_ARRAY_TASK_ID - 1 ))
CELL_TYPE="${CELL_TYPES[$IDX]}"

echo "============================================"
echo "Job:       $SLURM_JOB_ID (array task $SLURM_ARRAY_TASK_ID)"
echo "Cell type: $CELL_TYPE"
echo "Node:      $(hostname)"
echo "Started:   $(date)"
echo "============================================"

# ---------------------------------------------------------------------------
# Locate the .tbl file for this cell type
# ---------------------------------------------------------------------------
TBL_FILE=$(ls "$RESULTS_DIR"/"${CELL_TYPE}"_meta_analysis_*.tbl 2>/dev/null | \
           grep -v "[0-9]\.tbl$" | head -1)

if [ -z "$TBL_FILE" ]; then
  echo "ERROR: No .tbl file found for cell type $CELL_TYPE in $RESULTS_DIR"
  exit 1
fi
echo "TBL file: $TBL_FILE"

# ---------------------------------------------------------------------------
# Step 1: Extract loci (skip if already done)
# ---------------------------------------------------------------------------
LOCI_FILE="$COLOC_LOCI_DIR/${CELL_TYPE}_loci.tsv"
LOCUS_DATA="$COLOC_LOCI_DIR/${CELL_TYPE}_locus_data.tsv.gz"

if [ ! -f "$LOCI_FILE" ] || [ ! -f "$LOCUS_DATA" ]; then
  echo "--- Step 1: Extracting loci ---"
  Rscript "$SCRIPTS_DIR/01_extract_loci.R" \
    --cell_type   "$CELL_TYPE" \
    --tbl_file    "$TBL_FILE" \
    --output_dir  "$COLOC_LOCI_DIR" \
    --p_threshold 1e-5 \
    --window_kb   500
else
  echo "--- Step 1: Loci already extracted, skipping ---"
  NLOCI=$(tail -n +2 "$LOCI_FILE" | wc -l)
  echo "  Found $NLOCI loci in $LOCI_FILE"
fi

# ---------------------------------------------------------------------------
# Step 2: Run coloc (skip if already done)
# ---------------------------------------------------------------------------
COLOC_OUT="$COLOC_RES_DIR/${CELL_TYPE}_coloc_results.tsv"

if [ -f "$COLOC_OUT" ]; then
  echo "--- Step 2: Coloc results exist, skipping ($COLOC_OUT) ---"
  echo "  (Delete to rerun: rm $COLOC_OUT)"
else
  echo "--- Step 2: Running coloc.abf + coloc.susie ---"
  Rscript "$SCRIPTS_DIR/03_run_coloc.R" \
    --cell_type    "$CELL_TYPE" \
    --loci_file    "$LOCI_FILE" \
    --locus_data   "$LOCUS_DATA" \
    --disease_dir  "$DISEASE_DIR" \
    --output_dir   "$COLOC_RES_DIR" \
    --cell_type_N  "$CELL_TYPE_N" \
    --run_susie    TRUE \
    --pp_h4_susie  0.5 \
    --min_snps     50
fi

echo "============================================"
echo "Finished: $(date)"
echo "============================================"
