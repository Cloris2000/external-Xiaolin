#!/bin/bash
#SBATCH --job-name=fig6_v2_gws
#SBATCH --account=rrg-shreejoy
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=192
#SBATCH --time=00:30:00
#SBATCH --output=/scratch/zhoux156/logs/downstream_v2/fig6_v2_gws_%j.out
#SBATCH --error=/scratch/zhoux156/logs/downstream_v2/fig6_v2_gws_%j.err

# Figure 6 (v2), option A: loci chosen from the CTP-GWAS itself rather than from
# colocalisation, so Figure 6 does not repeat Figure 4's coloc-driven regional panels.
#   - cell types with >= 5 'Used in MGP' markers only (1-2-marker types L4 IT, L5 ET,
#     L6b, PAX6, L5/6 NP, L6 CT are marker-gene expression traits; their GWS hits are
#     cis-eQTLs of the marker)
#   - every genome-wide significant locus (one per 1 Mb; best cell type by lead P),
#     then the strongest suggestive locus of a cell type not yet shown
#   Dry run selected: VIP x TMEM106B, VIP x chr8:82.35 Mb, Oligodendrocyte x chr7:9.30 Mb.
# The previous coloc-driven Figure 6 outputs are kept under figure6_v2_coloc_driven/.
# Trillium allocates whole nodes; do not add --mem.

set -euo pipefail
source "${SITE_ENV:-/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/site_env.sh}"

V2META="${RESULTS_ROOT}/meta_analysis_15cohorts_hg19_v2"
DS="${RESULTS_ROOT}/downstream_v2"
FIG="${DS}/manuscript_figures"
KEEP="${FIG}/figure6_v2_coloc_driven"
echo "fig6_v2_gws  $(date -Is)"

mkdir -p "${KEEP}"
for f in figure6_cohort_heterogeneity_v2.png figure6_cohort_heterogeneity_v2.pdf figure6_v2_loci.tsv figure6_v2_per_cohort_effects.tsv; do
    if [ -e "${FIG}/${f}" ] && [ ! -e "${KEEP}/${f}" ]; then mv "${FIG}/${f}" "${KEEP}/${f}" && echo "  kept ${f} -> figure6_v2_coloc_driven/"; fi
done

cd "${FIG}"
"${RSCRIPT}" "${NF_DIR}/manuscript_figure/figure6_cohort_heterogeneity_v2.R" \
    --meta_dir "${V2META}" --metadata "${NF_DIR}/docs/meta_cohort_metadata.tsv" --out_dir "${FIG}" \
    --from_gws --loci_dir "${DS}/coloc/loci_full" --n_loci 3 --min_markers 5

echo ""; echo "Done  $(date -Is)"
cut -f1,2,6,7,8,10,11,18,19 "${FIG}/figure6_v2_loci.tsv" | column -t
