#!/bin/bash
#SBATCH --job-name=downstream_v2_redo2
#SBATCH --account=rrg-shreejoy
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=192
#SBATCH --time=03:00:00
#SBATCH --output=/scratch/zhoux156/logs/downstream_v2/downstream_v2_redo2_%j.out
#SBATCH --error=/scratch/zhoux156/logs/downstream_v2/downstream_v2_redo2_%j.err

# Second redo for the corrected 15-cohort meta (v2):
#
#   disease GWAS build   AD_Bellenguez2022, LBD_Chia2021 and PD_Nalls2019 "_hg19.tsv" files
#                        are actually GRCh38 (GWAS-Catalog harmonised coordinates copied
#                        through unlifted).  They are lifted to hg19 here into
#                        ${DS}/disease_gwas_hg19/; BD/MDD/SCZ (PGC daner, native hg19,
#                        rs429358 at its hg19 position) are linked unchanged.
#   coloc                re-run against the lifted set (six primary GWAS); previous v2
#                        outputs kept as coloc/*_grch38mix
#   figure 3             re-run with the Panel A lead tables now built in-script
#   figures 4 / coloc-rearranged / 6   depend on coloc -> re-run
#
# Trillium allocates whole nodes; do not add --mem.

set -uo pipefail
source "${SITE_ENV:-/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/site_env.sh}"

V2META="${RESULTS_ROOT}/meta_analysis_15cohorts_hg19_v2"
DS="${RESULTS_ROOT}/downstream_v2"
FIG="${DS}/manuscript_figures"
LOGD="${LOG_ROOT}/downstream_v2/${SLURM_JOB_ID:-manual}_redo2"
SCC_RES="${SCC_DIR}/results"
SRC_DIS="${SCC_RES}/coloc/disease_gwas"
DIS="${DS}/disease_gwas_hg19"
CHAIN="${REFS_ROOT}/reference_data/hg38ToHg19.over.chain.gz"
LIFT=(AD_Bellenguez2022 LBD_Chia2021 PD_Nalls2019)       # GRCh38 in disguise
KEEP=(BD_bip2024 MDD_MDD2025 SCZ_Trubetskoy2022)         # native hg19
DISEASES="MDD_MDD2025,BD_bip2024,SCZ_Trubetskoy2022,AD_Bellenguez2022,LBD_Chia2021,PD_Nalls2019"
mkdir -p "${LOGD}" "${FIG}" "${DIS}"
echo "downstream_v2_redo2  $(date -Is)"

step() { local name="$1"; shift; local t0; t0=$(date +%s)
    if "$@" > "${LOGD}/${name}.log" 2>&1; then printf "  OK    %-28s %4ss\n" "${name}" "$(( $(date +%s) - t0 ))"
    else printf "  FAIL  %-28s %4ss  (see %s)\n" "${name}" "$(( $(date +%s) - t0 ))" "${LOGD}/${name}.log"; echo "${name}" >> "${DS}/.failed_steps_redo2"; fi; }
rm -f "${DS}/.failed_steps_redo2"

# ---- Phase 0: hg19 disease set ------------------------------------------------
echo ""; echo "Phase 0  disease GWAS -> ${DIS}  $(date -Is)"
for d in "${KEEP[@]}"; do
    mkdir -p "${DIS}/${d}"; ln -sfn "${SRC_DIS}/${d}/${d}_hg19.tsv" "${DIS}/${d}/${d}_hg19.tsv"; echo "  linked ${d} (native hg19)"
done
for d in "${LIFT[@]}"; do
    ( step "liftover_${d}" "${PYTHON}" "${NF_DIR}/scripts/coloc/liftover_disease_gwas.py" \
          --input "${SRC_DIS}/${d}/${d}_hg19.tsv" --chain "${CHAIN}" --output "${DIS}/${d}/${d}_hg19.tsv" ) &
done
wait
for d in "${LIFT[@]}"; do [ -s "${DIS}/${d}/${d}_hg19.tsv" ] && cat "${LOGD}/liftover_${d}.log" | sed 's/^/    /'; done
# sanity: APOE rs429358 must now sit at hg19 19:45411941 in every file that carries it
for d in "${LIFT[@]}" "${KEEP[@]}"; do
    printf "    %-20s APOE@hg19 rows=%s  APOE@hg38 rows=%s\n" "$d" \
        "$(awk -F'\t' '$2==19 && $3>=45411900 && $3<=45412000' "${DIS}/${d}/${d}_hg19.tsv" | wc -l)" \
        "$(awk -F'\t' '$2==19 && $3>=44908600 && $3<=44908800' "${DIS}/${d}/${d}_hg19.tsv" | wc -l)"
done
n_ok=0; for d in "${LIFT[@]}"; do [ -s "${DIS}/${d}/${d}_hg19.tsv" ] && n_ok=$((n_ok+1)); done
[ "${n_ok}" -eq 3 ] || { echo "liftover failed for $((3-n_ok)) file(s); aborting before coloc"; exit 1; }

# ---- Phase 1: coloc on the hg19 set + Figure 3 --------------------------------
echo ""; echo "Phase 1  $(date -Is)"
for d in coloc_results_full full; do
    if [ -d "${DS}/coloc/${d}" ] && [ ! -d "${DS}/coloc/${d}_grch38mix" ]; then
        mv "${DS}/coloc/${d}" "${DS}/coloc/${d}_grch38mix" && echo "  kept previous coloc/${d} as coloc/${d}_grch38mix"; fi
done
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
  export MF_META_DIR="${V2META}" MF_CACHE_DIR="${FIG}/gwas_cache" MF_OUT_DIR="${FIG}"
  rm -f "${FIG}/figure3_data_checkpoint.rds"     # cache (with snp) from redo #1 is reused
  cd "${FIG}"
  step fig3_gwas_discovery    "${RSCRIPT}" "${NF_DIR}/manuscript_figure/figure3_gwas_discovery.R" \
  && step fig3_multict_manhattan "${RSCRIPT}" "${NF_DIR}/manuscript_figure/test_multict_gws_manhattan.R"
) &
wait

# ---- Phase 2: coloc-dependent figures -----------------------------------------
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
echo "figures:"; ls -la --time-style=+%H:%M "${FIG}"/*.png 2>/dev/null | awk '{print "  "$6"  "$7}' | sed 's#.*/##'
if [ -s "${DS}/.failed_steps_redo2" ]; then echo "FAILED steps:"; sed 's/^/  /' "${DS}/.failed_steps_redo2"; exit 1; fi
