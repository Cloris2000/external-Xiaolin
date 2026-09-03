#!/bin/bash
#SBATCH --job-name=anc_lead_fx
#SBATCH --partition=compute
#SBATCH --time=6:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=logs/meta_sensitivity/ancestry_lead_effects_%j.out
#SBATCH --error=logs/meta_sensitivity/ancestry_lead_effects_%j.err
#SBATCH --account=rrg-shreejoy

# Extract pooled GW leads + ancestry effects, then forest plots.
# Run after ancestry EUR/AFR metas finish (19 *.tbl each).

set -euo pipefail
source "${SITE_ENV:-/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/site_env.sh}"
cd "${NF_DIR}"
mkdir -p logs/meta_sensitivity

EUR_N=$(ls results/meta_sensitivity/ancestry_specific_EUR/*1.tbl 2>/dev/null | wc -l)
AFR_N=$(ls results/meta_sensitivity/ancestry_specific_AFR/*1.tbl 2>/dev/null | wc -l)
echo "EUR *1.tbl: ${EUR_N}  AFR *1.tbl: ${AFR_N}"
if [[ "${EUR_N}" -lt 19 || "${AFR_N}" -lt 19 ]]; then
  echo "WARNING: ancestry metas incomplete; proceeding with available cell types."
fi

RSCRIPT="${RSCRIPT:-Rscript}"

python3 scripts/extract_ancestry_lead_effects.py \
  --project-dir "${NF_DIR}" \
  --p-thresh 5e-8 \
  --clump-kb 500 \
  --max-leads-per-ct 6

"${RSCRIPT}" scripts/plot_ancestry_lead_forests.R \
  --effects results/meta_sensitivity/ancestry_lead_effects/ancestry_lead_effects.tsv \
  --output_dir results/meta_sensitivity/ancestry_lead_effects/forests \
  --max_leads_per_ct 6 \
  --combined TRUE

echo "Done. Outputs under results/meta_sensitivity/ancestry_lead_effects/"
