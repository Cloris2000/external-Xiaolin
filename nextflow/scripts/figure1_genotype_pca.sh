#!/bin/bash
#
# Figure 1 genotype panels: 1000 Genomes PCA projection for the 15 bulk cohorts.
#
# Produces the per-sample PC table behind Figure 1 panels A (superpopulation),
# B (cohort) and D (per-cohort ancestry composition).
#
# Method (standard reference-projection PCA):
#   1. Lift the 1000G hg38 reference DOWN to hg19, so it matches the 11 hg19
#      cohorts.  The project only ships hg38ToHg19.over.chain.gz and every meta
#      in this repo is hg19, so hg19 is the common space.  The four hg38 cohorts
#      (ROSMAP_array, GTEx_v10, AMP_AD_Rush, AMP_AD_Mayo - see
#      validation/genome_build_audit_bulk.txt) are lifted down the same way.
#   2. Intersect on chr:pos:REF:ALT variant IDs, keeping only variants where the
#      cohort and the reference agree on REF/ALT orientation (strand-unambiguous
#      SNPs only; A/T and C/G sites are dropped rather than guessed).
#   3. Fit PCs on the 1000G reference alone (--freq counts + --pca allele-wts),
#      then project every cohort sample into that reference space with
#      --score ... variance-standardize.  Cohort samples never influence the
#      axes, so PC1/PC2 are comparable across all 15 cohorts.
#
# Nothing here touches the GWAS inputs: it only reads the QC'd genotypes and
# writes to its own output directory.
#
# Usage (must run inside a SLURM allocation, not on a login node):
#   bash scripts/figure1_genotype_pca.sh
#
set -euo pipefail

source "${SITE_ENV:-/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/site_env.sh}"

PLINK="${REFS_ROOT}/tools/plink2"
CHAIN="${REFS_ROOT}/reference_data/hg38ToHg19.over.chain.gz"
KG_DIR="${DATA_ROOT}/refpanels/1000G_hg38"
# pyliftover lives in the `test` env (python_env does not have it).
PY_LIFT="${PY_LIFT:-${HOME}/miniforge3/envs/test/bin/python}"

# Per-cohort QC'd genotypes from the current GWAS run (cohort_gwas_v2).
GENO_ROOT="${GENO_ROOT:-${SCRATCH_ROOT}/results/cohort_gwas_v2/results}"
OUT_DIR="${OUT_DIR:-${SCRATCH_ROOT}/results/figure1_genotype_pca}"
TMP_DIR="${OUT_DIR}/tmp"

# Memory/threads for plink2 inside a whole-node allocation.
PLINK_MEM="${PLINK_MEM:-120000}"
PLINK_THREADS="${PLINK_THREADS:-16}"

# Cohort -> genotype file basename.  The three HBCC platform cohorts all carry
# the CMC_HBCC basename; everything else matches its cohort name.
COHORTS=(
    "ROSMAP:ROSMAP"                   "ROSMAP_array:ROSMAP_array"
    "Mayo:Mayo"                       "MSBB:MSBB"
    "CMC_MSSM:CMC_MSSM"               "CMC_PENN:CMC_PENN"
    "CMC_PITT:CMC_PITT"               "GTEx_v10:GTEx_v10"
    "NABEC:NABEC"                     "NIMH_HBCC_1M:CMC_HBCC"
    "NIMH_HBCC_h650:CMC_HBCC"         "NIMH_HBCC_Omni5M:CMC_HBCC"
    "GVEX:GVEX"                       "AMP_AD_Rush:AMP_AD_Rush"
    "AMP_AD_Mayo:AMP_AD_Mayo"
)

# Cohorts whose coordinates are hg38 and must be lifted to hg19.
# Source: validation/genome_build_audit_bulk.txt (landmarks hg38=5/5).
HG38_COHORTS=" ROSMAP_array GTEx_v10 AMP_AD_Rush AMP_AD_Mayo "

mkdir -p "${OUT_DIR}" "${TMP_DIR}"
echo "[figure1-pca] output: ${OUT_DIR}"

# ---------------------------------------------------------------------------
# Step 1: LD-pruned common-SNP backbone from the 1000G reference, lifted to hg19
# ---------------------------------------------------------------------------
if [[ ! -f "${TMP_DIR}/kg_hg19.pgen" ]]; then
    echo "[figure1-pca] 1/6 selecting common autosomal SNPs from 1000G"
    # Biallelic autosomal SNPs, MAF>=5%, non-missing; then LD-prune.
    "${PLINK}" \
        --pfile "${KG_DIR}/all_hg38" \
        --psam  "${KG_DIR}/hg38_corrected.psam" \
        --chr 1-22 \
        --snps-only just-acgt \
        --max-alleles 2 \
        --maf 0.05 \
        --geno 0.01 \
        --memory "${PLINK_MEM}" --threads "${PLINK_THREADS}" \
        --indep-pairwise 1000 50 0.2 \
        --out "${TMP_DIR}/kg_prune"

    echo "[figure1-pca]     pruned-in variants: $(wc -l < "${TMP_DIR}/kg_prune.prune.in")"

    # Materialize the pruned reference (still hg38 coordinates).
    "${PLINK}" \
        --pfile "${KG_DIR}/all_hg38" \
        --psam  "${KG_DIR}/hg38_corrected.psam" \
        --extract "${TMP_DIR}/kg_prune.prune.in" \
        --memory "${PLINK_MEM}" --threads "${PLINK_THREADS}" \
        --make-pgen --out "${TMP_DIR}/kg_hg38_pruned"

    # The reference stores IDs as "1:pos:REF:ALT"; every cohort VCF in this
    # project uses "chr1:pos:REF:ALT".  figure1_liftover_pvar.py emits the
    # chr-prefixed form, so the reference is normalized on the way through.

    echo "[figure1-pca] 2/6 lifting 1000G backbone hg38 -> hg19"
    "${PY_LIFT}" "$(dirname "$0")/figure1_liftover_pvar.py" \
        --pvar  "${TMP_DIR}/kg_hg38_pruned.pvar" \
        --chain "${CHAIN}" \
        --out-pvar    "${TMP_DIR}/kg_hg19.pvar.new" \
        --out-keep    "${TMP_DIR}/kg_hg19.keep" \
        --unmapped-log "${OUT_DIR}/kg_hg38_unmapped.tsv"

    # Drop variants that failed liftover, then swap in hg19 IDs/positions.
    "${PLINK}" \
        --pfile "${TMP_DIR}/kg_hg38_pruned" \
        --extract "${TMP_DIR}/kg_hg19.keep" \
        --memory "${PLINK_MEM}" --threads "${PLINK_THREADS}" \
        --make-pgen --out "${TMP_DIR}/kg_hg19_tmp"
    cp "${TMP_DIR}/kg_hg19.pvar.new" "${TMP_DIR}/kg_hg19_tmp.pvar"

    # Re-sort into hg19 coordinate order (liftover can reorder positions).
    "${PLINK}" \
        --pfile "${TMP_DIR}/kg_hg19_tmp" \
        --sort-vars \
        --memory "${PLINK_MEM}" --threads "${PLINK_THREADS}" \
        --make-pgen --out "${TMP_DIR}/kg_hg19"

    echo "[figure1-pca]     hg19 backbone variants: $(grep -vc '^#' "${TMP_DIR}/kg_hg19.pvar")"
else
    echo "[figure1-pca] 1-2/6 reusing existing hg19 1000G backbone"
fi

# ---------------------------------------------------------------------------
# Step 2: lift the four hg38 cohorts to hg19 and restrict every cohort to the
#         backbone, so all 15 live on one variant set
# ---------------------------------------------------------------------------
echo "[figure1-pca] 3/6 preparing per-cohort genotypes on the backbone"
KG_IDS="${TMP_DIR}/kg_backbone_ids.txt"
grep -v '^#' "${TMP_DIR}/kg_hg19.pvar" | cut -f3 > "${KG_IDS}"

for entry in "${COHORTS[@]}"; do
    cohort="${entry%%:*}"; base="${entry##*:}"
    src="${GENO_ROOT}/${cohort}/${base}.QC.final"
    out="${TMP_DIR}/coh_${cohort}"
    [[ -f "${out}.pgen" ]] && { echo "[figure1-pca]     ${cohort}: cached"; continue; }
    [[ -f "${src}.pgen" ]] || { echo "[figure1-pca]     ${cohort}: MISSING ${src}.pgen" >&2; exit 1; }

    if [[ "${HG38_COHORTS}" == *" ${cohort} "* ]]; then
        # hg38 cohort: lift its variant IDs to hg19 first, then intersect.
        "${PY_LIFT}" "$(dirname "$0")/figure1_liftover_pvar.py" \
            --pvar  "${src}.pvar" \
            --chain "${CHAIN}" \
            --out-pvar    "${TMP_DIR}/${cohort}.hg19.pvar" \
            --out-keep    "${TMP_DIR}/${cohort}.hg19.keep" \
            --unmapped-log "${OUT_DIR}/${cohort}_hg38_unmapped.tsv" \
            --restrict-to "${KG_IDS}"

        "${PLINK}" --pfile "${src}" \
            --extract "${TMP_DIR}/${cohort}.hg19.keep" \
            --memory "${PLINK_MEM}" --threads "${PLINK_THREADS}" \
            --make-pgen --out "${TMP_DIR}/${cohort}_lift"
        cp "${TMP_DIR}/${cohort}.hg19.pvar" "${TMP_DIR}/${cohort}_lift.pvar"
        "${PLINK}" --pfile "${TMP_DIR}/${cohort}_lift" \
            --extract "${KG_IDS}" --sort-vars \
            --memory "${PLINK_MEM}" --threads "${PLINK_THREADS}" \
            --make-pgen --out "${out}"
    else
        # hg19 cohort: no liftover, but IDs still have to be put in the same
        # canonical chr{N}:{pos}:{REF}:{ALT} form as the backbone before the
        # --extract can match anything.
        "${PLINK}" --pfile "${src}" \
            --memory "${PLINK_MEM}" --threads "${PLINK_THREADS}" \
            --make-just-pvar --out "${TMP_DIR}/${cohort}_ids"
        "${PY_LIFT}" "$(dirname "$0")/figure1_normalize_ids.py" \
            --pvar "${TMP_DIR}/${cohort}_ids.pvar" \
            --out  "${TMP_DIR}/${cohort}_ids.norm.pvar"

        # Rebuild the cohort with normalized IDs, then restrict to the backbone.
        "${PLINK}" --pfile "${src}" \
            --memory "${PLINK_MEM}" --threads "${PLINK_THREADS}" \
            --make-pgen --out "${TMP_DIR}/${cohort}_norm"
        cp "${TMP_DIR}/${cohort}_ids.norm.pvar" "${TMP_DIR}/${cohort}_norm.pvar"

        "${PLINK}" --pfile "${TMP_DIR}/${cohort}_norm" \
            --extract "${KG_IDS}" --sort-vars \
            --memory "${PLINK_MEM}" --threads "${PLINK_THREADS}" \
            --make-pgen --out "${out}"
        rm -f "${TMP_DIR}/${cohort}_norm".{pgen,pvar,psam} \
              "${TMP_DIR}/${cohort}_ids".pvar "${TMP_DIR}/${cohort}_ids.norm.pvar"
    fi
    echo "[figure1-pca]     ${cohort}: $(grep -vc '^#' "${out}.pvar") backbone variants, $(grep -vc '^#' "${out}.psam") samples"
done

# ---------------------------------------------------------------------------
# Step 3: common variant set across reference + all cohorts, orientation-checked
# ---------------------------------------------------------------------------
echo "[figure1-pca] 4/6 intersecting variants across reference and 15 cohorts"
"${PY_LIFT}" "$(dirname "$0")/figure1_common_variants.py" \
    --ref-pvar "${TMP_DIR}/kg_hg19.pvar" \
    --cohort-pvars "${TMP_DIR}"/coh_*.pvar \
    --out "${TMP_DIR}/common_ids.txt" \
    --report "${OUT_DIR}/variant_intersection_report.tsv"
echo "[figure1-pca]     shared variants: $(wc -l < "${TMP_DIR}/common_ids.txt")"

# ---------------------------------------------------------------------------
# Step 4: fit PCs on the 1000G reference only
# ---------------------------------------------------------------------------
echo "[figure1-pca] 5/6 fitting reference PCA (1000G only)"
# The 1000G panel is 3,202 samples but only 2,583 founders: it includes trio
# offspring.  Related samples distort a PCA (the axes start describing family
# structure), so the reference is restricted to founders with --keep-founders.
# That also resolves plink2's refusal to run "--freq counts" while nonfounders
# are present.  All five superpopulations keep 353-681 founders, which is ample.
#
# vcols is pinned to chrom,ref,alt so the .eigenvec.allele layout is fixed and
# known:  1=CHROM 2=ID 3=REF 4=ALT 5=A1 6..15=PC1..PC10.  (The default vcols
# includes `maybeprovref`, a column that appears only conditionally - pinning
# avoids the score column numbers below shifting out from under us.)
"${PLINK}" \
    --pfile "${TMP_DIR}/kg_hg19" \
    --extract "${TMP_DIR}/common_ids.txt" \
    --keep-founders \
    --freq counts \
    --pca 10 allele-wts vcols=chrom,ref,alt \
    --memory "${PLINK_MEM}" --threads "${PLINK_THREADS}" \
    --out "${OUT_DIR}/kg_ref_pca"

# Fail loudly if the layout is not what the --score column numbers assume.
_hdr=$(head -1 "${OUT_DIR}/kg_ref_pca.eigenvec.allele")
_c2=$(echo "${_hdr}" | cut -f2); _c5=$(echo "${_hdr}" | cut -f5); _c6=$(echo "${_hdr}" | cut -f6)
if [[ "${_c2}" != "ID" || "${_c5}" != "A1" || "${_c6}" != "PC1" ]]; then
    echo "[figure1-pca] ERROR: unexpected .eigenvec.allele layout: ${_hdr}" >&2
    echo "[figure1-pca]        expected col2=ID col5=A1 col6=PC1" >&2
    exit 1
fi
echo "[figure1-pca]     allele-weight layout verified (ID=2, A1=5, PC1=6)"

# ---------------------------------------------------------------------------
# Step 5: project every cohort into the reference PC space
# ---------------------------------------------------------------------------
echo "[figure1-pca] 6/6 projecting cohorts into reference PC space"
# Project the reference itself too, so panel A points and cohort points are
# produced by an identical scoring step (avoids the eigenvec/score scaling
# mismatch that makes projected samples look shrunk toward the origin).
#
# --read-freq makes `variance-standardize` use the 1000G allele frequencies for
# every target, instead of each cohort's own.  Without it each cohort would be
# standardized by its own variance and the resulting PCs would not be on a
# common scale.
for tgt in kg "${COHORTS[@]}"; do
    if [[ "${tgt}" == "kg" ]]; then
        # Founders only, matching the samples the PCs were fit on, so the
        # panel A/B backdrop is the reference panel itself and not its
        # trio offspring as well.
        name="1000G"; pfile="${TMP_DIR}/kg_hg19"; keep_arg="--keep-founders"
    else
        # Cohorts carry no pedigree, so every sample is a founder here; the
        # flag is deliberately not applied to them.
        name="${tgt%%:*}"; pfile="${TMP_DIR}/coh_${name}"; keep_arg=""
    fi
    "${PLINK}" \
        --pfile "${pfile}" \
        --extract "${TMP_DIR}/common_ids.txt" \
        ${keep_arg} \
        --read-freq "${OUT_DIR}/kg_ref_pca.acount" \
        --score "${OUT_DIR}/kg_ref_pca.eigenvec.allele" 2 5 header-read \
                no-mean-imputation variance-standardize \
        --score-col-nums 6-15 \
        --memory "${PLINK_MEM}" --threads "${PLINK_THREADS}" \
        --out "${TMP_DIR}/proj_${name}"
    echo "[figure1-pca]     projected ${name}"
done

# ---------------------------------------------------------------------------
# Step 6: assemble the tidy per-sample PC table
# ---------------------------------------------------------------------------
"${HOME}/miniforge3/envs/python_env/bin/python" "$(dirname "$0")/figure1_assemble_pca.py" \
    --proj-dir "${TMP_DIR}" \
    --kg-psam  "${KG_DIR}/hg38_corrected.psam" \
    --eigenval "${OUT_DIR}/kg_ref_pca.eigenval" \
    --ancestry-assignments "${NF_DIR}/docs/ancestry_specific/ancestry_sample_assignments.tsv" \
    --out-samples "${OUT_DIR}/figure1_pca_samples.tsv" \
    --out-summary "${OUT_DIR}/figure1_pca_summary.tsv"

echo "[figure1-pca] done."
echo "  per-sample PCs : ${OUT_DIR}/figure1_pca_samples.tsv"
echo "  cohort summary : ${OUT_DIR}/figure1_pca_summary.tsv"
