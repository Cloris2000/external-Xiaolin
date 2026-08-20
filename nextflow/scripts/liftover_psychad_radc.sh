#!/bin/bash
# ============================================================================
# Liftover PsychAD RADC (Mathys) snRNA cell-proportion GWAS sumstats hg38 -> hg19
# ============================================================================
# The RADC genotypes come from the imputed PsychAD combined.vcf.gz (GRCh38),
# so REGENIE step2 sumstats are in hg38.  Every other sn / bulk CTP cohort is
# hg19/b37, so these must be lifted before meta-analysis or sn-vs-bulk
# concordance comparison.
#
# This mirrors modules/liftover_sumstats.nf: it runs scripts/liftover_sumstats.py
# (pyliftover) on each per-cell-type .regenie.raw_p using the UCSC chain file.
#
# Prereqs (already satisfied):
#   - pyliftover installed in the `test` conda env
#   - reference_data/hg38ToHg19.over.chain.gz present
#   - RADC step2 sumstats produced in results/sn_psychad_radc_hodge/regenie_step2/
#
# Usage:
#   bash scripts/liftover_psychad_radc.sh
# ============================================================================
set -euo pipefail

PROJECT="/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow"
COHORT="PsychAD_RADC_snRNA"
IN_DIR="${PROJECT}/results/sn_psychad_radc_hodge/regenie_step2"
OUT_DIR="${PROJECT}/results/sn_psychad_radc_hodge/regenie_step2_hg19"
LOG_DIR="${OUT_DIR}/liftover_logs"
CHAIN="${PROJECT}/reference_data/hg38ToHg19.over.chain.gz"
SCRIPT="${PROJECT}/scripts/liftover_sumstats.py"

# Activate env with pyliftover
source /nethome/kcni/xzhou/.anaconda3/etc/profile.d/conda.sh
conda activate test

mkdir -p "${OUT_DIR}" "${LOG_DIR}"

[ -f "${CHAIN}" ]  || { echo "ERROR: chain file missing: ${CHAIN}"; exit 1; }
[ -f "${SCRIPT}" ] || { echo "ERROR: liftover script missing: ${SCRIPT}"; exit 1; }

shopt -s nullglob
raw_p_files=( "${IN_DIR}"/${COHORT}_*_step2.regenie.raw_p )
if [ ${#raw_p_files[@]} -eq 0 ]; then
    echo "ERROR: no ${COHORT}_*_step2.regenie.raw_p files in ${IN_DIR}"
    echo "       (has the RADC GWAS finished? check squeue / combined_gwas.log)"
    exit 1
fi

echo "Lifting ${#raw_p_files[@]} RADC sumstats files hg38 -> hg19 ..."
for f in "${raw_p_files[@]}"; do
    base=$(basename "$f")
    # PsychAD_RADC_snRNA_<cell_type>_step2.regenie.raw_p -> extract <cell_type>
    cell_type=${base#${COHORT}_}
    cell_type=${cell_type%_step2.regenie.raw_p}

    out="${OUT_DIR}/${COHORT}_${cell_type}.hg19_lifted.regenie.raw_p"
    unmapped="${LOG_DIR}/${COHORT}_${cell_type}_hg38_unmapped.tsv"

    python3 "${SCRIPT}" \
        --input-file   "$f" \
        --cohort       "${COHORT}" \
        --chain-file   "${CHAIN}" \
        --output-file  "${out}" \
        --unmapped-log "${unmapped}"
done

echo ""
echo "Done. Lifted files: ${OUT_DIR}/${COHORT}_<cell_type>.hg19_lifted.regenie.raw_p"
echo "Unmapped logs:      ${LOG_DIR}/"
