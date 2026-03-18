#!/bin/bash
# Run GWAS Results Comparison

echo "================================================="
echo "Comparing GWAS Results"
echo "================================================="
echo "Previous analysis: /external/rprshnas01/netdata_kcni/stlab/Xiaolin/WGS/ROSMAP_joint_wgs_step2"
echo "Nextflow pipeline: /external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/ROSMAP/regenie_step2"
echo "================================================="
echo ""

cd /external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow

# Activate conda environment (if needed)
# source ~/.anaconda3/etc/profile.d/conda.sh
# conda activate test

Rscript compare_gwas_results.R \
    --old_dir /external/rprshnas01/netdata_kcni/stlab/Xiaolin/WGS/ROSMAP_joint_wgs_step2 \
    --new_dir /external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/ROSMAP/regenie_step2 \
    --output_dir /external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/ROSMAP/comparison \
    --maf_threshold 0.05 \
    --cell_types "Astrocyte,Microglia,Oligodendrocyte,OPC,Endothelial,Pericyte,VLMC,IT,L4.IT,L5.ET,L5.6.IT.Car3,L5.6.NP,L6.CT,L6b,LAMP5,PAX6,PVALB,SST,VIP"

echo ""
echo "================================================="
echo "Comparison complete!"
echo "================================================="
echo "Results saved to: /external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/ROSMAP/comparison/"
echo ""
echo "Files generated:"
echo "  - *_manhattan_comparison.pdf : Side-by-side Manhattan plots"
echo "  - *_correlation.pdf : Scatter plots showing P-value correlation"
echo "  - comparison_statistics.csv : Summary statistics"
echo ""

