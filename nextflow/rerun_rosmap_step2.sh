#!/bin/bash
# Re-run Regenie Step 2 + convert_logp_to_p for 10 ROSMAP cell types
# whose Step1 pred.list was contaminated by sn_ROSMAP and has now been fixed.
# All 10 jobs use the April 7 bulk ROSMAP LOCO predictions.

set -euo pipefail

cd /external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow

REGENIE=/external/rprshnas01/kcni/mwainberg/software/regenie
PGEN=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/ROSMAP/ROSMAP.QC.final
PHENO=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/ROSMAP/phenotypes_RINT.txt
COVAR=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/ROSMAP/covariates.txt
STEP1=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/ROSMAP/regenie_step1
OUTDIR=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/ROSMAP/regenie_step2
LOGDIR=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/logs/ROSMAP_step2_rerun
THREADS=12
BSIZE=1000

mkdir -p "$LOGDIR"

CELL_TYPES=(Astrocyte Endothelial L4.IT L5.6.NP L5.ET L6.CT L6b OPC Oligodendrocyte VLMC)

for CT in "${CELL_TYPES[@]}"; do
    PRED="$STEP1/ROSMAP_${CT}_step1_pred.list"

    # Sanity check: pred.list must point to work/ROSMAP (not sn_ROSMAP)
    if ! grep -q "work/ROSMAP/" "$PRED"; then
        echo "ERROR: $PRED still points to sn_ROSMAP. Skipping $CT."
        continue
    fi

    # Remove old sentinel and output files so they are regenerated
    rm -f "$OUTDIR/ROSMAP_${CT}_step2_${CT}.regenie"
    rm -f "$OUTDIR/ROSMAP_${CT}_step2.regenie.raw_p"
    rm -f "$OUTDIR/ROSMAP_${CT}_step2.done"
    rm -f "$OUTDIR/ROSMAP_${CT}_logp_to_p.done"

    echo "Submitting Step2 + logp_to_p for: $CT"

    sbatch \
        --job-name="rosmap_s2_${CT}" \
        --partition=mediumtmp \
        --time=6:00:00 \
        --cpus-per-task=${THREADS} \
        --mem=48G \
        --output="${LOGDIR}/step2_${CT}_%j.out" \
        --error="${LOGDIR}/step2_${CT}_%j.err" \
        --wrap="
set -euo pipefail
cd ${OUTDIR}

echo '[Step2] Running Regenie Step 2 for ${CT}'
${REGENIE} \\
    --step 2 \\
    --threads ${THREADS} \\
    --verbose \\
    --pgen ${PGEN} \\
    --phenoFile ${PHENO} \\
    --covarFile ${COVAR} \\
    --bsize ${BSIZE} \\
    --pred ${PRED} \\
    --out ROSMAP_${CT}_step2

touch ROSMAP_${CT}_step2.done

echo '[logp_to_p] Converting LOG10P to P for ${CT}'
awk 'NR==1 {for(i=1;i<=NF;i++) if(\$i==\"LOG10P\") log10p_col=i; print \$0, \"P\"} NR>1 {print \$0, 10^(-\$log10p_col)}' \\
    ROSMAP_${CT}_step2_${CT}.regenie > ROSMAP_${CT}_step2.regenie.raw_p

touch ROSMAP_${CT}_logp_to_p.done

echo '[Done] ${CT}'
"
done

echo ""
echo "All 10 jobs submitted. Monitor with: squeue -u \$USER"
