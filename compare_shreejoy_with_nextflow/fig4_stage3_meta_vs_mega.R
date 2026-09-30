# Figure 4 -- Stage 3, combining cohorts.
#
# Mine runs REGENIE per cohort and combines with METAL, which yields an I2 per
# variant. His pools everyone into one mega-analysis with cohort as a covariate,
# which yields no I2 at all.
#
# The temptation is to read "meta gives 687 hits, mega gives 33" as a verdict on
# meta versus mega. It is not, and the control that shows why is his: he also
# ran a plain inverse-variance meta across his three ancestry strata, on the
# same people and the same variants. If pooling were doing the work, that meta
# would look like mine. It does not -- it recovers his two loci and nothing
# else. What separates the two pipelines at this stage is what went into them,
# not how the cohorts were combined.

source(file.path(
  "/project/rrg-shreejoy/zhoux156/external-Xiaolin/compare_shreejoy_with_nextflow",
  "theme_compare.R"))
suppressPackageStartupMessages(library(ggrepel))

disc <- read_stage("stage3_discovery.tsv")
mvm  <- read_stage("stage3_meta_vs_mega.tsv")
sub  <- read_stage("subclass_map.tsv")

# --- (a) my hits per trait, split by heterogeneity --------------------------
mine <- disc[pipeline == "mine"]
mine[, `:=`(n_gw = as.integer(n_gw), n_het = as.integer(n_gw_high_het))]
mine <- mine[n_gw > 0]

core_by_mgp <- sub[, .(n_core = sum(n_core)), by = mgp_trait]
mine <- merge(mine, core_by_mgp, by = "mgp_trait", all.x = TRUE)
mine[is.na(n_core), n_core := 0]

long <- rbind(
  mine[, .(trait, k = "high heterogeneity (I2 above threshold)", n = n_het)],
  mine[, .(trait, k = "low heterogeneity", n = n_gw - n_het)])
long[, trait := factor(trait, levels = mine[order(n_gw), trait])]

drop_lab <- mine[n_core == 0, trait]
pa <- ggplot(long, aes(n, trait, fill = k)) +
  geom_col(width = 0.72) +
  geom_text(data = mine, aes(x = n_gw, y = trait, label = n_gw),
            inherit.aes = FALSE, hjust = -0.25, size = 2.9,
            fontface = "bold", colour = "grey20") +
  scale_fill_manual(values = c("high heterogeneity (I2 above threshold)" =
                                 "#C1553B", "low heterogeneity" = "grey72"),
                    name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.16))) +
  scale_y_discrete(labels = function(x)
    ifelse(x %in% drop_lab, paste0(x, " *"), x)) +
  labs(x = "genome-wide significant variants (P < 5e-8)", y = NULL,
       title = "a  My meta-analysis: 687 hits, 529 of them high-I2",
       subtitle = paste0(
         "* marks a trait with no supertype that survived his portability ",
         "filter.\nThose nine traits carry 520 of the 687 hits.")) +
  theme_cmp(10) +
  theme(legend.position = "bottom")

# --- (b) his mega against his own ancestry meta -----------------------------
# Same people, same variants, two models. The only clean meta-vs-mega contrast
# available anywhere in either pipeline.
mega <- mvm[model == "mega (pooled)"]
meta <- mvm[model == "meta (ancestry IVW)"]
meta_loci <- fread(file.path(
  "/scratch/shreejoy/ctpgwas/results/meta_ancestry/loci.csv"))

loci_tbl <- data.table(
  model = factor(c("his mega\n(pooled, 5,008 people)",
                   "his meta\n(IVW over EUR/AFR/AMR)"),
                 levels = c("his mega\n(pooled, 5,008 people)",
                            "his meta\n(IVW over EUR/AFR/AMR)")),
  n_hits = c(nrow(mega), nrow(meta)),
  n_loci_gw = c(uniqueN(mega$locus), nrow(meta_loci)),
  n_loci_corrected = c(uniqueN(mega[as.numeric(P) < MEFF_LINE, locus]),
                       sum(as.numeric(meta_loci$best_P) < MEFF_LINE)),
  best_p = c(min(as.numeric(mega$P)), min(as.numeric(meta$P))))

# --- (b) heterogeneity of the hits, both pipelines --------------------------
# Not the same quantity on both sides, and the figure says so: mine is I2 over
# 15 cohorts each with its own phenotype pipeline, his is I2 over 3 ancestry
# strata of one harmonised phenotype. That difference is the point.
het_cmp <- data.table(
  who = factor(c("mine: 15 cohorts", "his: 3 ancestry strata"),
               levels = c("mine: 15 cohorts", "his: 3 ancestry strata")),
  pct = c(100 * mine[, sum(n_het) / sum(n_gw)],
          100 * mean(as.numeric(meta$I2) > 50)),
  n = c(mine[, sum(n_gw)], nrow(meta)))

pb <- ggplot(het_cmp, aes(who, pct, fill = who)) +
  geom_col(width = 0.5, alpha = 0.88) +
  geom_text(aes(label = sprintf("%.0f%%\nof %s hits", pct,
                                formatC(n, big.mark = ",", format = "d"))),
            vjust = -0.25, size = 3.1, fontface = "bold", lineheight = 1.1) +
  scale_fill_manual(values = c("mine: 15 cohorts" = PIPE_COL[["mine"]],
                               "his: 3 ancestry strata" = PIPE_COL[["his"]]),
                    guide = "none") +
  scale_y_continuous(limits = c(0, 100),
                     expand = expansion(mult = c(0, 0.08))) +
  labs(x = NULL, y = "% of significant hits with high heterogeneity",
       title = "b  Where the heterogeneity lives",
       subtitle = paste0(
         "Not the same quantity: mine is across 15 cohorts that each built ",
         "their own\nphenotype, his is across 3 ancestry strata of one ",
         "harmonised phenotype.\nThat is the difference, not an artefact of ",
         "the comparison.")) +
  theme_cmp(10)

# --- (c) the two models agree on the same loci ------------------------------
loci_long <- melt(
  loci_tbl, id.vars = "model",
  measure.vars = c("n_loci_gw", "n_loci_corrected"),
  variable.name = "threshold", value.name = "n_loci")
loci_long[, threshold := factor(
  threshold, levels = c("n_loci_gw", "n_loci_corrected"),
  labels = c("P < 5e-8", "P < 2.94e-09 (corrected)"))]

pc <- ggplot(loci_long, aes(model, n_loci, fill = threshold)) +
  geom_col(width = 0.6, position = position_dodge(0.66), alpha = 0.88) +
  geom_text(aes(label = n_loci), position = position_dodge(0.66),
            vjust = -0.35, size = 3.2, fontface = "bold") +
  scale_fill_manual(values = c("P < 5e-8" = "#7FA8BD",
                               "P < 2.94e-09 (corrected)" = PIPE_COL[["his"]]),
                    name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
  labs(x = NULL, y = "independent loci",
       title = "c  Meta and mega agree, on his data",
       subtitle = paste0(
         "At the corrected threshold both recover TMEM106B and GRN and ",
         "nothing else.\nPooling is not what separates the two pipelines at ",
         "this stage; the phenotype\nthat went into it is.")) +
  theme_cmp(10) +
  theme(legend.position = "bottom")

p <- (pa | (pb / pc)) + plot_layout(widths = c(1.05, 1)) +
  plot_annotation(
    title = "Stage 3: meta-analysis against mega-analysis",
    caption = paste(
      "Mine: results/meta_analysis_15cohorts_hg19_v2/heterogeneity/",
      "meta_heterogeneity_summary.tsv; high-I2 is that file's own threshold.",
      "\nHis mega: ctpgwas/results/loci/annotated_hits.csv, pooled arm.",
      "His meta: ctpgwas/results/meta_ancestry/{hits,loci}.csv,",
      "inverse-variance across EUR, AFR and AMR.",
      "\nPortability from celltype-composition/refs/taxonomy_DFC_2026.tsv."),
    theme = theme_cmp(12))

save_fig(p, "fig4_stage3_meta_vs_mega.png", width = 14.5, height = 9)

cat("\n")
print(loci_tbl)
cat(sprintf("\n  his ancestry-meta I2: median %.1f, mean %.1f, %% above 50 = %.1f\n",
            median(as.numeric(meta$I2)), mean(as.numeric(meta$I2)),
            100 * mean(as.numeric(meta$I2) > 50)))
