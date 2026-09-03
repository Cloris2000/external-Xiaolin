#!/usr/bin/env Rscript
# =============================================================================
# Extract per-donor genotype dosage + covariate-adjusted CTP for the two
# Panel C example variants, harmonized to the displayed effect (ALT) allele
# used in Panels A/B of the generalizability figure.
#
# For each contributing cohort:
#   * genotype: results/<cohort>/*.QC.final.{pgen,pvar,psam}  (plink2)
#       dosage counted for the effect (ALT) allele via --export-allele
#   * phenotype: results/<cohort>/phenotypes_RINT.txt  (RINT CTP = GWAS pheno)
#   * covariates: results/<cohort>/covariates.txt      (PC1-10, msex, age terms)
#   * join on IID (verified: overlap == genotype N in all cohorts)
#   * covariate-adjust: residualize RINT CTP on covariates (lm), then z-score
#     within cohort  ->  "Covariate-adjusted standardized CTP"
#
# No new association test is computed here (formal inference remains the
# cohort-level fixed-effect meta-analysis in Panels A/B/C forest).
# =============================================================================

suppressPackageStartupMessages({ library(data.table) })

ROOT   <- "/project/rrg-shreejoy/zhoux156/external-Xiaolin/nextflow"
GEN    <- file.path(ROOT, "results/meta_sensitivity/generalizability")
WORK   <- file.path(GEN, "genotype_ctp")
TAB    <- file.path(GEN, "figure_tables")
PLINK  <- "/project/rrg-shreejoy/pipeline_refs/tools/plink2"
dir.create(WORK, showWarnings = FALSE, recursive = TRUE)
dir.create(TAB,  showWarnings = FALSE, recursive = TRUE)

# cohort display label -> results subdir + pgen prefix (as used by the GWAS)
COHORTS <- list(
  list(lab = "ROSMAP (WGS)",  dir = "ROSMAP",           pre = "ROSMAP.QC.final"),
  list(lab = "Mayo",          dir = "Mayo",             pre = "Mayo.QC.final"),
  list(lab = "MSBB",          dir = "MSBB",             pre = "MSBB.QC.final"),
  list(lab = "CMC MSSM",      dir = "CMC_MSSM",         pre = "CMC_MSSM.QC.final"),
  list(lab = "CMC PENN",      dir = "CMC_PENN",         pre = "CMC_PENN.QC.final"),
  list(lab = "CMC PITT",      dir = "CMC_PITT",         pre = "CMC_PITT.QC.final"),
  list(lab = "NABEC",         dir = "NABEC",            pre = "NABEC.QC.final"),
  list(lab = "GVEX",          dir = "GVEX",             pre = "GVEX.QC.final"),
  list(lab = "HBCC (Omni5M)", dir = "NIMH_HBCC_Omni5M", pre = "CMC_HBCC.QC.final"),
  list(lab = "HBCC (1M)",     dir = "NIMH_HBCC_1M",     pre = "CMC_HBCC.QC.final"),
  list(lab = "HBCC (h650)",   dir = "NIMH_HBCC_h650",   pre = "CMC_HBCC.QC.final"))

# Panel C examples. The displayed EFFECT ALLELE is the one the meta pipeline
# harmonized to (NOT necessarily ALT): lead_effects_long "status" == "flipped"
# means the effect allele is REF; "ok" means the effect allele is ALT.
EXAMPLES <- list(
  list(cell = "Microglia", pheno = "Microglia", variant = "chr16:31298939:T:G",
       rsID = "rs9937837"),
  list(cell = "L5.ET",     pheno = "L5.ET",     variant = "chr2:10862188:G:A",
       rsID = "rs12468286"))

long_eff <- fread(file.path(GEN, "lead_effects_long.tsv"))
eff_allele_of <- function(vid) {
  st <- long_eff[lead_variant == vid & level == "cohort" & !is.na(beta), status]
  p <- strsplit(vid, ":", fixed = TRUE)[[1]]
  if (length(st) && mean(st == "flipped") > 0.5) p[3] else p[4]  # REF if flipped else ALT
}

geno_out  <- list()   # per-donor plot data
audit_out <- list()   # harmonization audit

for (ex in EXAMPLES) {
  vid <- ex$variant; eff <- eff_allele_of(vid)
  ref_alt <- strsplit(vid, ":", fixed = TRUE)[[1]]           # chr pos REF ALT
  REF <- ref_alt[3]; ALT <- ref_alt[4]
  # write per-example allele file (count the effect/ALT allele)
  af <- file.path(WORK, "allele.txt"); sf <- file.path(WORK, "snp.txt")
  writeLines(paste(vid, eff, sep = "\t"), af)
  writeLines(vid, sf)

  for (co in COHORTS) {
    pre <- file.path(ROOT, "results", co$dir, co$pre)
    if (!file.exists(paste0(pre, ".pgen"))) next
    stem <- file.path(WORK, sprintf("g_%s_%s", gsub("[^A-Za-z0-9]", "", co$lab),
                                    gsub("[^A-Za-z0-9]", "", vid)))
    rc <- system2(PLINK, c("--pfile", pre, "--extract", sf, "--export", "A",
                           "--export-allele", af, "--out", stem),
                  stdout = FALSE, stderr = FALSE)
    raw <- paste0(stem, ".raw")
    if (rc != 0 || !file.exists(raw)) {
      audit_out[[length(audit_out) + 1]] <- data.table(
        cell_type = ex$cell, variant = vid, displayed_effect_allele = eff,
        original_genotype_coding = paste0("REF=", REF, ";ALT=", ALT),
        harmonized_dosage_coding = paste0("copies of ", eff),
        cohort = co$lab, donors_before_qc = 0L, donors_after_qc = 0L,
        n_0 = 0L, n_1 = 0L, n_2 = 0L,
        exclusions_reasons = "variant absent in cohort genotypes")
      next
    }
    g <- fread(raw)
    dcol <- grep(paste0(gsub("([][{}().^$*+?|\\\\])", "\\\\\\1", vid), "_"),
                 names(g), value = TRUE)
    if (!length(dcol)) dcol <- names(g)[ncol(g)]
    setnames(g, dcol, "dosage")
    g <- g[, .(IID, dosage)]
    n_before <- nrow(g)

    ph <- fread(file.path(ROOT, "results", co$dir, "phenotypes_RINT.txt"))
    if (!(ex$pheno %in% names(ph))) next
    ph <- ph[, c("IID", ex$pheno), with = FALSE]; setnames(ph, ex$pheno, "ctp")
    cv <- fread(file.path(ROOT, "results", co$dir, "covariates.txt"))
    cvcols <- setdiff(names(cv), c("FID", "IID"))

    d <- merge(merge(g, ph, by = "IID"), cv[, c("IID", cvcols), with = FALSE], by = "IID")
    n_geno_pheno <- nrow(d)
    # exclude missing genotype / phenotype / covariates
    keep <- stats::complete.cases(d[, c("dosage", "ctp", cvcols), with = FALSE])
    excl_missing <- sum(!keep)
    d <- d[keep]

    # covariate-adjust (drop constant covariate columns) then z-score within cohort
    use_cv <- cvcols[sapply(cvcols, function(c) length(unique(d[[c]])) > 1)]
    fml <- as.formula(paste("ctp ~", paste(use_cv, collapse = " + ")))
    d[, adj := residuals(lm(fml, data = d))]
    d[, adj_z := (adj - mean(adj)) / sd(adj)]
    d[, dgrp := pmin(pmax(round(dosage), 0L), 2L)]

    geno_out[[length(geno_out) + 1]] <- data.table(
      cell_type = ex$cell, rsID = ex$rsID, variant = vid,
      effect_allele = eff, cohort = co$lab, IID = d$IID,
      dosage = d$dosage, dosage_group = d$dgrp, adj_ctp_z = d$adj_z)

    cnt <- tabulate(d$dgrp + 1L, nbins = 3)
    audit_out[[length(audit_out) + 1]] <- data.table(
      cell_type = ex$cell, variant = vid, displayed_effect_allele = eff,
      original_genotype_coding = paste0("REF=", REF, ";ALT=", ALT,
                                        "; plink counted ", eff),
      harmonized_dosage_coding = paste0("copies of ", eff, " (0/1/2)"),
      cohort = co$lab, donors_before_qc = n_geno_pheno,
      donors_after_qc = nrow(d), n_0 = cnt[1], n_1 = cnt[2], n_2 = cnt[3],
      exclusions_reasons = if (excl_missing > 0)
        sprintf("%d dropped (missing genotype/phenotype/covariate)", excl_missing)
      else "none")
  }
}

geno  <- rbindlist(geno_out,  use.names = TRUE)
audit <- rbindlist(audit_out, use.names = TRUE)
fwrite(geno,  file.path(TAB, "panelC_genotype_CTP_data.tsv"), sep = "\t")
fwrite(audit, file.path(TAB, "panelC_genotype_harmonization_audit.tsv"), sep = "\t")

# pooled genotype-group counts + <10 safeguard
pooled <- geno[, .(n = .N), by = .(cell_type, variant, dosage_group)]
setorder(pooled, cell_type, variant, dosage_group)
cat("Per-example pooled genotype-group counts:\n"); print(pooled)
flag <- pooled[n < 10]
if (nrow(flag)) { cat("\n*** WARNING: genotype groups with <10 pooled donors ***\n"); print(flag) }
cat("\nDONE extract_genotype_ctp\n")
