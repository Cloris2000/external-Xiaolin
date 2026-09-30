#!/bin/bash
#SBATCH --job-name=downstream_v2_redo3
#SBATCH --account=rrg-shreejoy
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=192
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/zhoux156/logs/downstream_v2/downstream_v2_redo3_%j.out
#SBATCH --error=/scratch/zhoux156/logs/downstream_v2/downstream_v2_redo3_%j.err

# Third redo for the corrected 15-cohort meta (v2): coloc + its dependent figures only.
#
#   03_run_coloc.R / 05_regional_coloc_plots.R matched cell-type and disease SNPs on
#   the full chr:pos:ref:alt ID, so a disease file whose effect allele is written
#   in the opposite order (PD_Nalls2019: 0 exact vs 2,375 swapped matches at a test
#   locus) never merged and PD had 23 tests vs ~460 for the others.  Matching is now
#   on an allele-order-agnostic key with the existing exact-or-flipped allele check.
#   Inputs: v2 loci (unchanged) and the hg19-lifted disease set from redo #2.
#   Previous coloc outputs kept as coloc/*_idmerge.  Also re-runs the multi-cell-type
#   Manhattan (svglite guard).
# Trillium allocates whole nodes; do not add --mem.

set -uo pipefail
source "${SITE_ENV:-/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/site_env.sh}"

V2META="${RESULTS_ROOT}/meta_analysis_15cohorts_hg19_v2"
DS="${RESULTS_ROOT}/downstream_v2"
FIG="${DS}/manuscript_figures"
DIS="${DS}/disease_gwas_hg19"
LOGD="${LOG_ROOT}/downstream_v2/${SLURM_JOB_ID:-manual}_redo3"
DISEASES="MDD_MDD2025,BD_bip2024,SCZ_Trubetskoy2022,AD_Bellenguez2022,LBD_Chia2021,PD_Nalls2019"
mkdir -p "${LOGD}"
echo "downstream_v2_redo3  $(date -Is)"
[ -s "${DIS}/PD_Nalls2019/PD_Nalls2019_hg19.tsv" ] || { echo "lifted disease set missing: ${DIS}"; exit 1; }

step() { local name="$1"; shift; local t0; t0=$(date +%s)
    if "$@" > "${LOGD}/${name}.log" 2>&1; then printf "  OK    %-28s %4ss\n" "${name}" "$(( $(date +%s) - t0 ))"
    else printf "  FAIL  %-28s %4ss  (see %s)\n" "${name}" "$(( $(date +%s) - t0 ))" "${LOGD}/${name}.log"; echo "${name}" >> "${DS}/.failed_steps_redo3"; fi; }
rm -f "${DS}/.failed_steps_redo3"

for d in coloc_results_full full; do
    if [ -d "${DS}/coloc/${d}" ] && [ ! -d "${DS}/coloc/${d}_idmerge" ]; then
        mv "${DS}/coloc/${d}" "${DS}/coloc/${d}_idmerge" && echo "  kept previous coloc/${d} as coloc/${d}_idmerge"; fi
done

echo ""; echo "Phase 1  $(date -Is)"
(
  export META_DIR="${V2META}" COLOC_LOCI_DIR="${DS}/coloc/loci_full" COLOC_RES_DIR="${DS}/coloc/coloc_results_full"
  export DISEASE_DIR="${DIS}" CELL_TYPE_N=3554 COLOC_DISEASES="${DISEASES}"
  step coloc_all_loci bash "${NF_DIR}/run_coloc_all_loci.sbatch" \
  && step coloc_summarize "${RSCRIPT}" "${NF_DIR}/scripts/coloc/04_summarize_coloc.R" \
        --results_dir "${DS}/coloc/coloc_results_full" --loci_dir "${DS}/coloc/loci_full" --output_dir "${DS}/coloc/full" \
  && step coloc_regional_plots "${RSCRIPT}" "${NF_DIR}/scripts/coloc/05_regional_coloc_plots.R" \
        --sig_file "${DS}/coloc/full/coloc_significant.tsv" --loci_dir "${DS}/coloc/loci_full" \
        --disease_dir "${DIS}" --meta_dir "${V2META}" --output_dir "${DS}/coloc/full/plots/regional" --genome_build 37
) &
(
  export MF_OUT_DIR="${FIG}"; cd "${FIG}"
  step fig3_multict_manhattan "${RSCRIPT}" "${NF_DIR}/manuscript_figure/test_multict_gws_manhattan.R"
) &
wait

echo ""; echo "Phase 2  $(date -Is)"
MF="${NF_DIR}/manuscript_figure"
export MF_OUT_DIR="${FIG}" MF_META_DIR="${V2META}" MF_AUTO_REGIONAL=1
export MF_COLOC_DIR="${DS}/coloc/full" MF_LOCI_DIR="${DS}/coloc/loci_full" MF_DISEASE_DIR="${DIS}"
cd "${FIG}"
step fig4_coloc_validation  "${RSCRIPT}" "${MF}/figure4_coloc_validation.R"
step fig_coloc_rearranged   "${RSCRIPT}" "${MF}/figure_coloc_rearranged.R"
step fig6_heterogeneity_v2  "${RSCRIPT}" "${MF}/figure6_cohort_heterogeneity_v2.R" \
    --meta_dir "${V2META}" --metadata "${NF_DIR}/docs/meta_cohort_metadata.tsv" --out_dir "${FIG}" \
    --from_coloc "${DS}/coloc/full/coloc_all_results.tsv" --loci_dir "${DS}/coloc/loci_full" --n_loci 3

echo ""; echo "Done  $(date -Is)"
echo "coloc tests per disease:"; awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next}{d[$h["disease"]]++} END{for(k in d) printf "  %-22s %d\n",k,d[k]}' "${DS}/coloc/full/coloc_all_results.tsv" 2>/dev/null | sort
echo "coloc events (PP.H4>=0.5): $(( $(wc -l < "${DS}/coloc/full/coloc_significant.tsv" 2>/dev/null || echo 1) - 1 ))"
if [ -s "${DS}/.failed_steps_redo3" ]; then echo "FAILED steps:"; sed 's/^/  /' "${DS}/.failed_steps_redo3"; exit 1; fi
