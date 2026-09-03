#!/bin/bash
#SBATCH --job-name=rosmap_manhattan_comp
#SBATCH --output=/scratch/zhoux156/logs/comparison_%j.out
#SBATCH --error=/scratch/zhoux156/logs/comparison_%j.err
#SBATCH --time=2:00:00
#SBATCH --partition=compute
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --account=rrg-shreejoy

echo "================================================="
echo "Comparing GWAS Results (SLURM Job)"
echo "Job ID: $SLURM_JOB_ID"
echo "Start time: $(date)"
echo "================================================="
echo "Previous analysis: /external/rprshnas01/netdata_kcni/stlab/Xiaolin/WGS/ROSMAP_joint_wgs_step2"
echo "Nextflow pipeline: /project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/results/ROSMAP/regenie_step2"
echo "Output directory: /project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/results/ROSMAP/comparison"
echo "================================================="
echo ""

cd /project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow

# Create logs directory if it doesn't exist
mkdir -p /project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/logs

# Activate conda environment if needed
# source ~/miniforge3/etc/profile.d/conda.sh
# conda activate test

# Run comparison
Rscript compare_gwas_results.R \
    --old_dir /external/rprshnas01/netdata_kcni/stlab/Xiaolin/WGS/ROSMAP_joint_wgs_step2 \
    --new_dir /project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/results/ROSMAP/regenie_step2 \
    --output_dir /project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/results/ROSMAP/comparison \
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
echo "Results saved to: /project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/results/ROSMAP/comparison/"
echo ""
echo "Files generated:"
echo "  - *_manhattan_comparison.pdf : Side-by-side Manhattan plots"
echo "  - *_correlation.pdf : Scatter plots showing P-value correlation"
echo "  - comparison_statistics.csv : Summary statistics"
echo ""

exit $EXIT_CODE

