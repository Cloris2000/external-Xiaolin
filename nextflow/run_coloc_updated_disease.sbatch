#!/bin/bash
#SBATCH --job-name=coloc_bd_mdd_upd
#SBATCH --array=1-19%5
#SBATCH --time=4:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/coloc/coloc_upd_%A_%a.out
#SBATCH --error=logs/coloc/coloc_upd_%A_%a.err

# Re-run coloc.abf for updated BD (bip2024) and MDD (MDD2025) only,
# then merge into coloc_results_full replacing Mullins/Wray rows.
#
# Prerequisites:
#   bash scripts/coloc/02_download_disease_gwas.sh  (standardizes BD_bip2024, MDD_MDD2025)
#   results/coloc/loci_full/ from prior full coloc run
#
# After all array jobs finish:
#   Rscript scripts/coloc/04_summarize_coloc.R \
#     --results_dir results/coloc/coloc_results_full \
#     --loci_dir    results/coloc/loci_full \
#     --output_dir  results/coloc/full

set -euo pipefail

export R_LIBS="/external/rprshnas01/kcni/xzhou/.anaconda3/lib/R/library"
export PATH="/nethome/kcni/xzhou/.anaconda3/bin:$PATH"

PROJECT_DIR="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow"
COLOC_LOCI_DIR="$PROJECT_DIR/results/coloc/loci_full"
COLOC_EXISTING="$PROJECT_DIR/results/coloc/coloc_results_full"
COLOC_NEW="$PROJECT_DIR/results/coloc/coloc_results_updated"
DISEASE_DIR="$PROJECT_DIR/results/coloc/disease_gwas"
SCRIPTS_DIR="$PROJECT_DIR/scripts/coloc"
CELL_TYPE_N=2812
MIN_SNPS=10
DISEASE_INCLUDE="BD_bip2024,MDD_MDD2025"
REPLACE_DISEASES="BD_Mullins2021,MDD_Wray2018"

mkdir -p logs/coloc "$COLOC_NEW"

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
echo "Diseases:  $DISEASE_INCLUDE"
echo "Started:   $(date)"
echo "============================================"

LOCI_FILE="$COLOC_LOCI_DIR/${CELL_TYPE}_loci.tsv"
LOCUS_DATA="$COLOC_LOCI_DIR/${CELL_TYPE}_locus_data.tsv.gz"

if [ ! -f "$LOCI_FILE" ] || [ ! -f "$LOCUS_DATA" ]; then
  echo "ERROR: Missing loci files for $CELL_TYPE in $COLOC_LOCI_DIR"
  exit 1
fi

for f in "$DISEASE_DIR/BD_bip2024/BD_bip2024_hg19.tsv" \
         "$DISEASE_DIR/MDD_MDD2025/MDD_MDD2025_hg19.tsv"; do
  if [ ! -s "$f" ]; then
    echo "ERROR: Missing standardized disease file: $f"
    echo "Run: bash scripts/coloc/02_download_disease_gwas.sh"
    exit 1
  fi
done

echo "--- coloc.abf for updated BD/MDD ---"
Rscript "$SCRIPTS_DIR/03_run_coloc.R" \
  --cell_type         "$CELL_TYPE" \
  --loci_file         "$LOCI_FILE" \
  --locus_data        "$LOCUS_DATA" \
  --disease_dir       "$DISEASE_DIR" \
  --output_dir        "$COLOC_NEW" \
  --cell_type_N       "$CELL_TYPE_N" \
  --run_susie         FALSE \
  --min_snps          "$MIN_SNPS" \
  --disease_include   "$DISEASE_INCLUDE"

echo "--- merge into coloc_results_full ---"
Rscript "$SCRIPTS_DIR/merge_coloc_results.R" \
  --existing_dir      "$COLOC_EXISTING" \
  --new_file          "$COLOC_NEW/${CELL_TYPE}_coloc_results.tsv" \
  --replace_diseases  "$REPLACE_DISEASES" \
  --output_dir        "$COLOC_EXISTING"

echo "============================================"
echo "Finished: $(date)"
echo "============================================"
