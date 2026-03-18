#!/bin/bash
#SBATCH --job-name=rosmap_manhattan_comp
#SBATCH --output=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/logs/comparison_%j.out
#SBATCH --error=/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/logs/comparison_%j.err
#SBATCH --time=2:00:00
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

echo "================================================="
echo "Comparing GWAS Results (SLURM Job)"
echo "Job ID: $SLURM_JOB_ID"
echo "Start time: $(date)"
echo "================================================="
echo "Previous analysis: /external/rprshnas01/netdata_kcni/stlab/Xiaolin/WGS/ROSMAP_joint_wgs_step2"
echo "Nextflow pipeline: /external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/ROSMAP/regenie_step2"
echo "Output directory: /external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/ROSMAP/comparison"
echo "================================================="
echo ""

cd /external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow

# Create logs directory if it doesn't exist
mkdir -p /external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/logs

# Activate conda environment if needed
# source ~/.anaconda3/etc/profile.d/conda.sh
# conda activate test

# Run comparison
Rscript compare_gwas_results.R \
    --old_dir /external/rprshnas01/netdata_kcni/stlab/Xiaolin/WGS/ROSMAP_joint_wgs_step2 \
    --new_dir /external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/ROSMAP/regenie_step2 \
    --output_dir /external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/ROSMAP/comparison \
    --maf_threshold 0.05 \
    --cell_types "Astrocyte,Microglia,Oligodendrocyte,OPC,Endothelial,Pericyte,VLMC,IT,L4.IT,L5.ET,L5.6.IT.Car3,L5.6.NP,L6.CT,L6b,LAMP5,PAX6,PVALB,SST,VIP"

EXIT_CODE=$?

echo ""
echo "================================================="
if [ $EXIT_CODE -eq 0 ]; then
    echo "Comparison complete successfully!"
else
    echo "Comparison failed with exit code: $EXIT_CODE"
fi
echo "End time: $(date)"
echo "================================================="
echo "Results saved to: /external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/ROSMAP/comparison/"
echo ""
echo "Files generated:"
echo "  - *_manhattan_comparison.pdf : Side-by-side Manhattan plots"
echo "  - *_correlation.pdf : Scatter plots showing P-value correlation"
echo "  - comparison_statistics.csv : Summary statistics"
echo ""

exit $EXIT_CODE

