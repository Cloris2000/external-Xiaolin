# Figure 5 -- Stage 4, colocalisation with disease.
#
# This stage is close to a controlled experiment. His refs/disease_gwas/panel.tsv
# points AD_Bellenguez2022, PD_Nalls2019 and LBD_Chia2021 at
# /scratch/zhoux156/results/downstream_v2/disease_gwas_hg19/, annotated
# "Xiaolin's standardised hg19". For those three the disease side is literally
# the same files, so any difference in the posterior comes from the cell-type
# side alone.
#
# The two designs differ in shape as well as size. Mine tests every cell type
# whose scan produced a locus, 12 of them at TMEM106B. His tests one cell type
# per locus, the one that discovered it. So the comparison is his single point
# against my spread, not a heatmap.

source(file.path(
  "/project/rrg-shreejoy/zhoux156/external-Xiaolin/compare_shreejoy_with_nextflow",
  "theme_compare.R"))
suppressPackageStartupMessages(library(ggrepel))

col <- read_stage("stage4_coloc.tsv")
col[, PP.H4 := as.numeric(PP.H4)]

# Diseases both panels can show. AD, PD and LBD are the same underlying files.
SHARED <- c("AD", "PD", "LBD", "SCZ", "BD", "MDD")
SAME_FILE <- c("AD", "PD", "LBD")
col[, disease := factor(disease, levels = SHARED)]

locus_panel <- function(loc, title, subtitle) {
  d <- col[locus == loc & !is.na(disease)]
  mine <- d[pipeline == "mine"]
  his  <- d[pipeline == "his"]

  p <- ggplot(mapping = aes(disease, PP.H4)) +
    geom_hline(yintercept = 0.8, linetype = "22", colour = "grey40",
               linewidth = 0.4) +
    scale_x_discrete(drop = FALSE, limits = SHARED,
                     labels = function(x)
                       ifelse(x %in% SAME_FILE, paste0(x, "\u2020"), x)) +
    scale_y_continuous(limits = c(0, 1.08), breaks = seq(0, 1, 0.25)) +
    labs(x = NULL, y = "PP.H4 (one shared causal variant)",
         title = title, subtitle = subtitle) +
    theme_cmp(10)

  if (nrow(mine)) {
    p <- p +
      geom_jitter(data = mine, width = 0.13, height = 0, size = 1.9,
                  colour = PIPE_COL[["mine"]], alpha = 0.75) +
      geom_text_repel(
        data = mine[, .SD[which.max(PP.H4)], by = disease][PP.H4 > 0.8],
        aes(label = sprintf("%s %.2f", cell_type, PP.H4)),
        size = 2.5, colour = PIPE_COL[["mine"]], fontface = "bold",
        direction = "y", nudge_y = 0.07, segment.colour = "grey75")
  } else {
    p <- p + annotate(
      "text", x = 3.5, y = 0.5, size = 3.2, colour = PIPE_COL[["mine"]],
      fontface = "bold", lineheight = 1.2,
      label = paste("my pipeline never tested this locus:",
                    "\nno cell type reached significance there,",
                    "\nso it was never passed to coloc"))
  }

  p + geom_point(data = his, shape = 23, size = 3.4, stroke = 1.1,
                 fill = PIPE_COL[["his"]], colour = "black") +
    geom_text_repel(data = his[PP.H4 > 0.5],
                    aes(label = sprintf("%s %.3f", cell_type, PP.H4)),
                    size = 2.5, colour = PIPE_COL[["his"]], fontface = "bold",
                    direction = "y", nudge_y = -0.12,
                    segment.colour = "grey75")
}

pa <- locus_panel(
  "TMEM106B", "a  TMEM106B: where both pipelines test, both agree",
  paste0("AD 0.915 against 0.902, PD 0.033 against 0.041, LBD 0.066 against ",
         "0.070, SCZ 0.023\nagainst 0.021. The coloc step itself is not where ",
         "the pipelines diverge. MDD is the one\ngap, and the two panels do ",
         "not read the same MDD file."))

pb <- locus_panel(
  "GRN/FAM171A2", "b  GRN / FAM171A2",
  paste0("His strongest colocalisation anywhere, PP.H4 0.9999 with AD on ",
         "Sst_19,\nfine-mapped to a single variant at posterior 0.990."))

# --- (c) how much was tested, and how much survived -------------------------
scale_tbl <- col[, .(n_tests = .N, n_strong = sum(PP.H4 > 0.8, na.rm = TRUE)),
                 by = pipeline]
scale_tbl[, who := factor(c(mine = "mine", his = "his")[pipeline],
                          levels = c("mine", "his"))]
scale_long <- melt(scale_tbl, id.vars = "who",
                   measure.vars = c("n_tests", "n_strong"),
                   variable.name = "k", value.name = "n")
scale_long[, k := factor(k, levels = c("n_tests", "n_strong"),
                         labels = c("coloc tests run",
                                    "with PP.H4 > 0.8"))]

pc <- ggplot(scale_long, aes(who, n, fill = who)) +
  geom_col(width = 0.52, alpha = 0.88) +
  geom_text(aes(label = formatC(n, big.mark = ",", format = "d")),
            vjust = -0.3, size = 3.1, fontface = "bold") +
  facet_wrap(~ k, scales = "free_y") +
  scale_fill_manual(values = PIPE_COL, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(x = NULL, y = NULL,
       title = "c  Scale of the colocalisation stage",
       subtitle = paste0(
         "Mine tests every cell type at every locus its scan produced, so the ",
         "count inherits\nstage 3's 687 hits. His tests the discovering cell ",
         "type at 9 tiered loci.")) +
  theme_cmp(10)

p <- ((pa | pb) / pc) + plot_layout(heights = c(1.3, 1)) +
  plot_annotation(
    title = "Stage 4: colocalisation, with the disease side held fixed",
    subtitle = paste("Red points are my cell types, blue diamonds his.",
                     "Dagger marks a disease where both pipelines read the",
                     "same summary-statistics file."),
    caption = paste(
      "Mine: /scratch/zhoux156/results/downstream_v2/coloc/coloc_results_full/",
      "*_coloc_results.tsv, coloc.abf.",
      "\nHis: ctpgwas/results/coloc/coloc_results.tsv; panel from",
      "ctp-gwas/refs/disease_gwas/panel.tsv, which sources AD, PD and LBD from",
      "Xiaolin's standardised hg19 files.",
      "\nHis panel also carries ADHD, cognitive performance and IQ, which mine",
      "does not; those columns are omitted rather than shown as absent."),
    theme = theme_cmp(12))

save_fig(p, "fig5_stage4_coloc.png", width = 14, height = 10)

cat("\n")
print(scale_tbl[, .(pipeline, n_tests, n_strong)])
cat("\n  strongest per locus and pipeline:\n")
print(col[locus %in% c("TMEM106B", "GRN/FAM171A2") & disease %in% SHARED,
          .SD[which.max(PP.H4)], by = .(pipeline, locus, disease)][
            order(locus, disease, pipeline),
            .(pipeline, locus, disease, cell_type, PP.H4 = round(PP.H4, 3))])
