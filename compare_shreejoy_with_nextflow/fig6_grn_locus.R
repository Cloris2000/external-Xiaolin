# Figure 6 -- the GRN locus, the one his pipeline found and mine did not.
#
# GRN is the cleanest case in the comparison. It is his second established
# locus, it fine-maps to a single variant at posterior 0.990, and it
# colocalises with Alzheimer's at PP.H4 0.9999 using Xiaolin's own AD summary
# statistics. In the 19-class meta-analysis it is not significant in any cell
# type, so it never entered fine-mapping or colocalisation.
#
# The figure asks the narrower, answerable question: is the signal absent from
# my data, or present but diluted? Panel (c) answers it by looking at the
# effect estimate rather than the p value, and the answer is diluted -- my SST
# carries the same sign at roughly a quarter the magnitude, with the same
# standard error.

source(file.path(
  "/project/rrg-shreejoy/zhoux156/external-Xiaolin/compare_shreejoy_with_nextflow",
  "theme_compare.R"))
suppressPackageStartupMessages(library(ggrepel))

CHR <- 17
XLIM <- c(44.15e6, 44.60e6)           # hg38
GENE <- c(44345086, 44353106)         # GRN, hg38
HIS_LEAD <- "chr17:44352876:C:T"

mine <- read_stage("my_loci_hg38.tsv")[CHROM == CHR]
his  <- read_stage("his_loci.tsv")[CHROM == CHR & arm == "pooled"]
tmap <- read_stage("trait_map.tsv")

mine[, `:=`(P = as.numeric(P), POS = as.numeric(POS_hg38),
            BETA = as.numeric(BETA_ALT), SE = as.numeric(SE))]
his[, `:=`(P = as.numeric(P), POS = as.numeric(POS_hg38),
           BETA = as.numeric(BETA), SE = as.numeric(SE))]
mine <- mine[POS %between% XLIM & is.finite(P)]
his  <- his[POS %between% XLIM & is.finite(P)]

tmap_u <- unique(tmap[, .(trait = regenie_trait, class_label, mgp_trait)])
his <- merge(his, tmap_u, by = "trait", all.x = TRUE)
his[trait == "NeuronGliaRatio", `:=`(class_label = "NA", mgp_trait = "NA")]
his[, cls := fifelse(grepl("GABA", class_label), "Inhibitory",
             fifelse(grepl("Glutam", class_label), "Excitatory",
             fifelse(grepl("Non-neuronal", class_label), "Glia", "NA")))]

gene_band <- annotate("rect", xmin = GENE[1], xmax = GENE[2],
                      ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.6)
gene_lab <- annotate("text", x = mean(GENE), y = Inf, label = "GRN",
                     vjust = 1.4, size = 3, fontface = "italic",
                     colour = "grey30")
lead_line <- geom_vline(xintercept = 44352876, linetype = "22",
                        colour = "#B03A2E", linewidth = 0.4)

ymax <- max(-log10(c(mine$P, his$P)), na.rm = TRUE) * 1.08

# --- (a) mine ----------------------------------------------------------------
my_best <- mine[which.min(P)]
pa <- ggplot(mine, aes(POS, -log10(P))) +
  gene_band + gene_lab + lead_line +
  geom_point(colour = "grey70", size = 0.6, alpha = 0.55) +
  geom_point(data = mine[trait == "SST"], colour = PIPE_COL[["mine"]],
             size = 1.1) +
  geom_hline(yintercept = -log10(GW_LINE), linetype = "22",
             colour = "grey35", linewidth = 0.4) +
  geom_text_repel(data = my_best,
                  aes(label = sprintf("best of all 19 classes: %s, P = %.2g",
                                      trait, P)),
                  size = 2.9, fontface = "bold", colour = "grey20",
                  nudge_y = 1.2) +
  scale_x_continuous(labels = function(x) sprintf("%.2f Mb", x / 1e6)) +
  coord_cartesian(ylim = c(0, ymax)) +
  labs(x = NULL, y = expression(-log[10](P)),
       title = "a  My meta-analysis, all 19 classes (SST highlighted)",
       subtitle = paste("Nothing approaches significance, so this locus was",
                        "never passed to fine-mapping or coloc")) +
  theme_cmp(10)

# --- (b) his -----------------------------------------------------------------
pb <- ggplot(his, aes(POS, -log10(P))) +
  gene_band + gene_lab + lead_line +
  geom_point(aes(colour = cls), size = 0.6, alpha = 0.5) +
  geom_hline(yintercept = -log10(GW_LINE), linetype = "22",
             colour = "grey35", linewidth = 0.4) +
  geom_hline(yintercept = -log10(MEFF_LINE), colour = "grey20",
             linewidth = 0.4) +
  geom_text_repel(data = his[which.min(P)],
                  aes(label = sprintf("%s, P = %.2g", trait, P)),
                  size = 2.9, fontface = "bold", colour = "grey20",
                  nudge_y = 1.2) +
  scale_colour_manual(values = CLASS_COL, name = NULL) +
  scale_x_continuous(labels = function(x) sprintf("%.2f Mb", x / 1e6)) +
  coord_cartesian(ylim = c(0, ymax)) +
  labs(x = "chr17 position (hg38)", y = expression(-log[10](P)),
       title = "b  His mega-analysis, 33 supertypes",
       subtitle = paste("Red line marks his lead", HIS_LEAD,
                        "- fine-mapped to this variant at posterior 0.990")) +
  theme_cmp(10) + theme(legend.position = "bottom")

# --- (c) absent, or present but diluted? ------------------------------------
lead_row <- his[variant == HIS_LEAD][1]
stopifnot(!is.na(lead_row$POS))

his_lead <- his[variant == HIS_LEAD, .(label = trait, cls, mgp_trait,
                                       beta = BETA, se = SE, p = P,
                                       pipeline = "his")]
my_lead <- mine[POS == lead_row$POS & toupper(REF) == toupper(lead_row$REF) &
                toupper(ALT) == toupper(lead_row$ALT),
                .(label = trait, cls = "NA", mgp_trait = trait,
                  beta = BETA, se = SE, p = P, pipeline = "mine")]

both <- rbind(his_lead, my_lead)
both[, `:=`(lo = beta - 1.96 * se, hi = beta + 1.96 * se)]
both[, grp := ifelse(pipeline == "mine", "my 19 broad classes",
                     "his 33 supertypes")]
both[, lab_f := factor(label, levels = both[order(grp, beta), label])]

pc <- ggplot(both, aes(beta, lab_f, colour = pipeline)) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.45,
                 alpha = 0.7) +
  geom_point(size = 1.9) +
  geom_point(data = both[p < 5e-8], size = 3.1, shape = 21, stroke = 0.9,
             fill = NA, colour = "black") +
  facet_grid(grp ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_colour_manual(values = PIPE_COL, guide = "none") +
  labs(x = sprintf("effect at %s, with 95%% CI", HIS_LEAD), y = NULL,
       title = "c  Diluted, not absent",
       subtitle = paste0(
         "Same allele on both sides; ringed points clear 5e-8. My SST points ",
         "the same way as\nhis Sst supertypes, at -0.033 against their mean ",
         "-0.093 and his best -0.136.\nThe standard errors are the same ",
         "(0.022 against 0.021), so what is lost is effect\nsize, not power. ",
         "A four-fold attenuation is enough to cost the locus entirely.")) +
  theme_cmp(10) +
  theme(axis.text.y = element_text(size = rel(0.72)),
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 90),
        panel.grid.major.y = element_line(colour = "grey94"))

p <- ((pa / pb) | pc) + plot_layout(widths = c(1.2, 1)) +
  plot_annotation(
    title = "GRN: the locus his pipeline found and mine never saw",
    caption = paste(
      "Mine: results/meta_analysis_15cohorts_hg19_v2/*.annotated.tsv, hg19",
      "positions carried to hg38 by +1,922,632 (data/coord_offsets.tsv).",
      "\nHis: ctpgwas/results/step2/all_shrunk, pooled arm.",
      "Fine-mapping posterior and coloc from",
      "ctp-gwas/stages/report/gwas_report/REPORT.md sections 4 and 5."),
    theme = theme_cmp(12))

save_fig(p, "fig6_grn_locus.png", width = 15, height = 9.5)

cat("\n  my best in window:", my_best$trait, "P =", signif(my_best$P, 3), "\n")
cat("  his best in window:", his[which.min(P), trait], "P =",
    signif(his[, min(P)], 3), "\n\n")
out <- both[order(pipeline, p),
            .(pipeline, label, beta = round(beta, 4), se = round(se, 4),
              p = signif(p, 3))]
print(head(out, 14))
fwrite(out, file.path(DATA_DIR, "fig6_grn_lead_effects.tsv"), sep = "\t")
