#!/bin/bash
# Compare one cohort's cohort_gwas_v2 output against the SCC run.
#
#   compare_cohort_to_scc.sh <cohort> <new_dir> <scc_dir> <out_dir>
#
# Writes into <out_dir>:
#   gwas_all_ct.tsv        per cell type: n_new n_old n_shared r_beta r_se r_logp lead SNPs   (compare_raw_p.R)
#   deconv.tsv             cell_proportions / phenotypes_RINT vs SCC                          (compare_deconv_to_scc.R)
#   orientation.tsv        copy of the deconv orientation gate result
#   pca.tsv                |r| of PC1..PC10 between new and SCC pca.csv on shared samples
#   samples.tsv            QC'd sample counts new / old / shared (psam)
#   variants_per_chr.tsv   QC'd variant counts per chromosome new vs old (pvar) with % change
#   flags.txt              anomalies: r_beta<0.9, shared<95%, lead SNP change, sample or
#                          per-chr variant count change >5%, negative orientation
# Exit status is 0 even when flags are raised; the flags are the report.
#
# Expected differences (documented, not anomalies): REGENIE step-1 SNP filter is
# hom_alt here vs mac1 in the SCC run; deconvolution uses exclude_bio_from_tech_cov;
# NABEC chr1 is the recovered file (SCC used the truncated one).

set -uo pipefail
cohort="$1"; new="$2"; old="$3"; out="$4"
source "${SITE_ENV:-/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow/site_env.sh}"
mkdir -p "$out"
FLAGS="$out/flags.txt"; : > "$FLAGS"
flag() { echo "$*" | tee -a "$FLAGS"; }

echo "== [$cohort] GWAS per cell type =="
rows=()
for f in "$new"/regenie_step2/*_step2.regenie.raw_p; do
    [ -s "$f" ] || continue
    b=$(basename "$f" _step2.regenie.raw_p); ct=${b#${cohort}_}
    o="$old/regenie_step2/$(basename "$f")"
    if [ ! -s "$o" ]; then flag "GWAS ${ct}: no SCC counterpart"; continue; fi
    "${RSCRIPT}" "${NF_DIR}/scripts/compare_raw_p.R" --cohort "$cohort" --cell_type "$ct" \
        --new "$f" --old "$o" --out "$out/.gwas_${ct}.tsv" >/dev/null 2>&1 || flag "GWAS ${ct}: compare_raw_p.R failed"
    rows+=("$out/.gwas_${ct}.tsv")
done
if [ ${#rows[@]} -gt 0 ]; then
    { head -1 "${rows[0]}"; for r in "${rows[@]}"; do [ -s "$r" ] && tail -n +2 "$r"; done; } > "$out/gwas_all_ct.tsv"
    rm -f "$out"/.gwas_*.tsv
    awk -F'\t' 'NR>1 {
        if ($6+0 < 0.9)            printf "GWAS %s: r_beta=%.3f < 0.9\n", $2, $6;
        if ($5/($4>0?$4:1) < 0.95) printf "GWAS %s: shared variants %.1f%% of SCC\n", $2, 100*$5/$4;
        if ($11 != "TRUE")         printf "GWAS %s: lead SNP changed %s -> %s\n", $2, $10, $9;
    }' "$out/gwas_all_ct.tsv" | tee -a "$FLAGS"
    echo "  $(($(wc -l < "$out/gwas_all_ct.tsv")-1)) cell types compared"
fi

echo "== [$cohort] deconvolution / phenotypes =="
"${RSCRIPT}" "${NF_DIR}/scripts/compare_deconv_to_scc.R" --cohort "$cohort" \
    --new_prop "$new/cell_proportions.csv" --old_prop "$old/cell_proportions.csv" \
    --new_pheno "$new/phenotypes_RINT.txt" --old_pheno "$old/phenotypes_RINT.txt" \
    --out "$out/deconv.tsv" >/dev/null 2>&1 || flag "DECONV: compare_deconv_to_scc.R failed"
[ -s "$new/deconv_orientation.tsv" ] && cp "$new/deconv_orientation.tsv" "$out/orientation.tsv" && \
    awk -F'\t' 'NR>1 && $5=="NEGATIVE"{printf "ORIENTATION %s: NEGATIVE marker correlation r=%s (n_markers=%s)\n",$2,$4,$3}
                NR>1 && $5=="WEAK"    {printf "ORIENTATION %s: weak marker correlation r=%s (n_markers=%s) - estimate unreliable in this cohort\n",$2,$4,$3}' \
        "$out/orientation.tsv" | tee -a "$FLAGS"

echo "== [$cohort] PCA =="
np="$new/pca.csv"; op="$old/pca.csv"
if [ -s "$np" ] && [ -s "$op" ]; then
"${RSCRIPT}" -e '
suppressMessages(library(data.table))
a <- fread("'"$np"'"); b <- fread("'"$op"'"); setnames(a,1,"id"); setnames(b,1,"id"); a[,id:=as.character(id)]; b[,id:=as.character(id)]
m <- merge(a, b, by="id", suffixes=c(".new",".old"))
pcs <- intersect(grep("^PC",names(a),value=TRUE), grep("^PC",names(b),value=TRUE))[1:10]; pcs <- pcs[!is.na(pcs)]
res <- rbindlist(lapply(pcs, function(p) data.table(cohort="'"$cohort"'", pc=p, n_shared=nrow(m), abs_r=round(abs(cor(m[[paste0(p,".new")]], m[[paste0(p,".old")]])),4))))
fwrite(res, "'"$out/pca.tsv"'", sep="\t"); print(res)
' 2>&1 | tail -12
awk -F'\t' 'NR>1 && $2 ~ /^PC[1-3]$/ && $4 < 0.9 {printf "PCA %s: |r|=%s vs SCC\n",$2,$4}' "$out/pca.tsv" | tee -a "$FLAGS"
else echo "  pca.csv missing (new: $([ -s "$np" ] && echo yes || echo no), old: $([ -s "$op" ] && echo yes || echo no))"; fi

echo "== [$cohort] QC'd samples and variants =="
npsam=$(ls "$new"/*.QC.final.psam 2>/dev/null | head -1); opsam=$(ls "$old"/*.QC.final.psam 2>/dev/null | head -1)
if [ -s "$npsam" ] && [ -s "$opsam" ]; then
    n1=$(awk 'NR>1' "$npsam" | wc -l); n2=$(awk 'NR>1' "$opsam" | wc -l)
    ns=$(comm -12 <(awk 'NR>1{print $2}' "$npsam" | sort) <(awk 'NR>1{print $2}' "$opsam" | sort) | wc -l)
    printf "cohort\tn_new\tn_old\tn_shared\n%s\t%s\t%s\t%s\n" "$cohort" "$n1" "$n2" "$ns" > "$out/samples.tsv"
    echo "  samples new=$n1 old=$n2 shared=$ns"
    awk -v a="$n1" -v b="$n2" 'BEGIN{ if (b>0 && (a-b)/b > 0.05 || (b-a)/b > 0.05) printf "SAMPLES: %d new vs %d SCC\n", a, b }' | tee -a "$FLAGS"
fi
npvar=$(ls "$new"/*.QC.final.pvar 2>/dev/null | head -1); opvar=$(ls "$old"/*.QC.final.pvar 2>/dev/null | head -1)
if [ -s "$npvar" ] && [ -s "$opvar" ]; then
    awk -v C="$cohort" 'BEGIN{OFS="\t"} FNR==1{fi++} !/^#/{ if(fi==1) a[$1]++; else b[$1]++ }
        END{ print "cohort","chr","n_new","n_old","pct_change";
             for(c=1;c<=22;c++){ pc=(b[c]>0)?100*(a[c]-b[c])/b[c]:"NA"; print C,c,a[c]+0,b[c]+0,(pc=="NA"?pc:sprintf("%.1f",pc)) } }' \
        "$npvar" "$opvar" > "$out/variants_per_chr.tsv"
    awk -F'\t' 'NR>1 && $5!="NA" && ($5+0 > 5 || $5+0 < -5) {printf "VARIANTS chr%s: %s new vs %s SCC (%s%%)\n",$2,$3,$4,$5}' "$out/variants_per_chr.tsv" | tee -a "$FLAGS"
    echo "  variants: $(awk -F'\t' 'NR>1{a+=$3;b+=$4} END{printf "new=%d old=%d", a, b}' "$out/variants_per_chr.tsv")"
fi

echo "== [$cohort] flags: $(wc -l < "$FLAGS") =="
exit 0
