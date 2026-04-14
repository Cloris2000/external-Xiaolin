# Multi-Ancestry GWAS Meta-Analysis: Checklist vs Current Pipeline

This document compares the standard multi-ancestry meta-analysis workflow (as used in high-impact studies and consortia) with the current Nextflow pipeline (`combined_pipeline_v2.nf`) and METAL meta-analysis (`standalone_meta_analysis.nf` + `modules/metal_meta_analysis.nf`).

---

## 1. What Your Pipeline Already Does Well ✅

| Requirement | Your pipeline | Status |
|-------------|----------------|--------|
| **Meta (not mega)** | Per-cohort GWAS → summary stats → METAL meta | ✅ Appropriate |
| **GWAS per cohort** | REGENIE step1 + step2 per cohort | ✅ |
| **Same phenotype definition** | Same cell-type proportions (RINT) across cohorts | ✅ |
| **Covariates** | Age, sex, top PCs, (batch where applicable) | ✅ |
| **Tool** | REGENIE | ✅ Common (with BOLT/PLINK/SAIGE) |
| **Meta scheme** | METAL `SCHEME STDERR` (inverse-variance weighted) | ✅ Preferred in modern studies |
| **Population stratification** | Top PCs included (all PCs from plink2 `--pca`, typically 20) | ✅ (10–20 is standard; you use all from PCA) |
| **MAF** | `maf_threshold` in genotype QC (e.g. 0.05 for HBCC; 0.01 in base config) | ✅ Stricter than “MAF > 0.01” |
| **Imputation QC (where used)** | `mach_r2_filter` in genotyping_qc (e.g. 0.8 for ROSMAP, CMC, MSBB, Mayo) | ✅ When set in config |
| **Post-meta QC** | Manhattan + QQ + lambda (genomic inflation) in `generate_meta_manhattan.R` | ✅ |
| **Allele frequency tracking** | METAL `AVERAGEFREQ ON`, `MINMAXFREQ ON` | ✅ Helps spot AF inconsistencies |

---

## 2. Discrepancies and Gaps

### 2.1 Harmonization (Before Meta) — **Missing**

Standard practice is a dedicated harmonization step before METAL:

| Item | Standard practice | Your pipeline | Discrepancy |
|------|-------------------|---------------|-------------|
| **Strand-ambiguous SNPs** | Remove or fix A/T and C/G when strand is unclear | No step; METAL gets REGENIE output as-is | **No explicit removal of strand-ambiguous SNPs.** METAL can align strands for non–strand-ambiguous SNPs; for A/T and C/G, METAL expects consistent strand across cohorts. Without harmonization, ambiguous SNPs can cause wrong effect directions. |
| **Effect-allele alignment** | Align effect allele (and REF/ALT) across cohorts | METAL aligns using ALLELE0/ALLELE1; no pre-check | **No pre-METAL check** that effect alleles and allele coding are consistent (e.g. same genome build, same REF/ALT). Relies on METAL’s internal alignment. |
| **Allele frequency consistency** | Check AF across cohorts; flag or drop large discrepancies | METAL reports MinFreq/MaxFreq only; no filter | **No filter** on AF discrepancy (e.g. drop if MaxFreq − MinFreq > threshold). You only inspect Min/Max in METAL output. |
| **INFO (imputation quality)** | Filter INFO &lt; 0.8 (or similar) in sumstats or in QC | No INFO column in REGENIE output; filter only at genotype QC via `mach_r2_filter` | **Cohort-dependent:** Where `mach_r2_filter` is set (e.g. 0.8 for ROSMAP, CMC, Mayo, MSBB), INFO-style QC is applied at genotype level. **HBCC configs do not set `mach_r2_filter`**, so for HBCC no explicit INFO/R² filter is applied. |
| **Rare / unstable variants** | Often filter very rare (e.g. MAF &lt; 0.01) or require present in ≥2 cohorts | MAF at genotype QC only; no “present in ≥2 cohorts” | **No “present in ≥2 cohorts”** rule in meta. Variants in a single cohort are still meta-analyzed. |
| **Genome build** | Same build (e.g. GRCh37/hg19) and consistent variant IDs | Not enforced in pipeline; assumed same build per cohort | **No explicit check** that all cohorts use the same build and matching variant IDs. |

**Recommendation:** Add a harmonization step (e.g. MungeSumstats, EasyQC, or a custom script) between per-cohort `.regenie.raw_p` and METAL: same build, remove strand-ambiguous (or fix strand), optional INFO filter if you add INFO to sumstats, and optionally “present in ≥2 cohorts” and AF-consistency filters.

---

### 2.2 METAL Meta-Analysis — **Partial**

| Item | Standard practice | Your pipeline | Discrepancy |
|------|-------------------|---------------|-------------|
| **Heterogeneity** | Report Cochran’s Q and I² per variant | METAL script uses `ANALYZE` only | **Heterogeneity not requested.** METAL supports `ANALYZE HETEROGENEITY` to output Q and I². Your `metal_meta_analysis.nf` does not use it, so heterogeneity stats are not produced. |
| **Weighting** | Inverse-variance (STDERR) or sample-size | `SCHEME STDERR` | ✅ Correct |
| **FLIP** | Often ON to allow METAL to flip strands | `FLIP OFF` | **Strand alignment** is effectively “no flipping” in METAL; alignment relies on allele names. With a proper harmonization step (same build, unambiguous alleles), FLIP OFF is fine; without harmonization, you depend entirely on METAL’s allele matching. |

**Recommendation:** In `metal_meta_analysis.nf`, add `ANALYZE HETEROGENEITY` (and keep `ANALYZE` if you want both fixed-effect and heterogeneity in one run, per METAL docs) so that meta output includes Q and I² for reporting and for sensitivity (e.g. random-effects or ancestry-specific reporting when I² is high).

---

### 2.3 GWAS Model and Covariates — **Mostly Aligned**

| Item | Standard practice | Your pipeline | Discrepancy |
|------|-------------------|---------------|-------------|
| **Same covariates across cohorts** | Same set (e.g. age, sex, 10–20 PCs, batch) | Per-cohort `covariates.txt` (PCA + clinical + optional batch); structure is cohort-specific | **Not formally enforced.** Covariates are consistent by design (same scripts and config), but there is no automated check that column names and number of PCs are identical across cohorts. |
| **Number of PCs** | Often 10–20 | All PCs from plink2 `--pca` (typically 20) | ✅ Within typical range |

---

### 2.4 Post-Meta QC and Reporting — **Partial**

| Item | Standard practice | Your pipeline | Discrepancy |
|------|-------------------|---------------|-------------|
| **Manhattan plot** | Yes | `generate_meta_manhattan.R` | ✅ |
| **QQ plot** | Yes | `generate_meta_manhattan.R` | ✅ |
| **Genomic inflation (lambda)** | Yes | Lambda computed in `generate_meta_manhattan.R` | ✅ |
| **LD score regression intercept** | Often reported to guard against residual stratification | Not implemented | **No LDSC intercept.** Only lambda is computed. For publication-grade multi-ancestry meta, LDSC intercept is often requested by reviewers. |

**Recommendation:** Optionally add an LDSC step on meta sumstats (e.g. `ldsc --h2` with appropriate LD scores for your build) and report the intercept.

---

### 2.5 Optional / Advanced (Not Required by Checklist)

- **Ancestry-aware meta (MANTRA, MR-MEGA):** Not in your pipeline; checklist says these are optional; standard fixed-effect IVW meta is fine.
- **Random-effects meta:** Not implemented; consider only if you routinely see high I² and want a sensitivity analysis.
- **Sample overlap correction:** METAL supports `OVERLAP ON` for SCHEME SAMPLESIZE; you use SCHEME STDERR, so this is less relevant. If some individuals appear in multiple cohorts, consider documenting or adjusting (e.g. cohort design) rather than changing scheme.

---

## 3. Summary Table: Fulfilled vs Gaps

| Category | Fulfilled | Gaps |
|----------|-----------|------|
| **Pipeline structure** | Per-cohort GWAS → meta; REGENIE; METAL STDERR | — |
| **Harmonization** | — | No pre-METAL harmonization; no strand-ambiguous removal; no INFO in sumstats; HBCC no mach_r2; no “≥2 cohorts” or AF-consistency filter; no formal genome-build check |
| **METAL** | STDERR, AVERAGEFREQ, MINMAXFREQ | No `ANALYZE HETEROGENEITY` (no Q / I²) |
| **Covariates** | Same design (age, sex, PCs, batch) | No automated cross-cohort consistency check |
| **Post-meta** | Manhattan, QQ, lambda | No LDSC intercept |

---

## 4. Recommended Fixes (Priority)

1. **High (for multi-ancestry reporting):**  
   - Add **`ANALYZE HETEROGENEITY`** in METAL so meta results include Q and I².

2. **High (for correctness):**  
   - Add a **harmonization step** before METAL: same build, remove (or fix) strand-ambiguous SNPs, optional INFO filter if you add INFO to sumstats.  
   - For **HBCC**, set **`mach_r2_filter`** (e.g. 0.8) in the relevant configs if imputed data has an R2/INFO field, so low-quality imputed variants are dropped at genotype QC.

3. **Medium:**  
   - Optionally filter meta results to variants **present in ≥2 cohorts** (or document that single-cohort variants are included).  
   - Optionally add **LDSC intercept** on meta sumstats for reporting.

4. **Low:**  
   - Optional automated check that **covariate columns** (and PC count) are consistent across cohorts.

---

## 5. References in Your Repo

- **GWAS:** `gwas_pipeline.nf`, `modules/regenie_step1.nf`, `modules/regenie_step2.nf`
- **Log10P → P:** `modules/convert_logp_to_p.nf` (produces `.regenie.raw_p`)
- **METAL:** `modules/metal_meta_analysis.nf` (METAL script: SCHEME STDERR, no HETEROGENEITY)
- **Standalone meta:** `standalone_meta_analysis.nf` (reads existing `.regenie.raw_p`, runs METAL)
- **Post-meta plots:** `scripts/generate_meta_manhattan.R` (Manhattan, QQ, lambda)
- **Genotype QC (MAF, HWE, R2):** `scripts/genotyping_qc.py`; `mach_r2_filter` in `modules/genotyping_qc_step.nf` and cohort configs
