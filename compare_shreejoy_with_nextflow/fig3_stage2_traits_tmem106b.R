# Figure 3 -- Stage 2, phenotype resolution, at TMEM106B.
#
# The two pipelines scan the same haplotype with different phenotypes. Mine
# averages a whole class; his keeps supertypes apart. The report says the
# inhibitory class is not unanimous there -- 13 supertypes move one way and 3
# the other -- so a class average cancels. This figure tests that claim against
# both sets of summary statistics rather than restating it.
#
# Coordinates: my hg19 positions are carried to hg38 by the offset derived in
# 04_align_coords.py, so both pipelines share one axis.

source(file.path(
  "/project/rrg-shreejoy/zhoux156/external-Xiaolin/compare_shreejoy_with_nextflow",
  "theme_compare.R"))
suppressPackageStartupMessages(library(ggrepel))

CHR <- 7
XLIM <- c(12.00e6, 12.45e6)           # hg38, centred on TMEM106B
GENE <- c(12211240, 12243367)         # TMEM106B, hg38

mine <- read_stage("my_loci_hg38.tsv")[CHROM == CHR]
his  <- read_stage("his_loci.tsv")[CHROM == CHR & arm == "pooled"]
tmap <- read_stage("trait_map.tsv")

mine[, `:=`(P = as.numeric(P), POS = as.numeric(POS_hg38),
            BETA = as.numeric(BETA_ALT), SE = as.numeric(SE))]
his[, `:=`(P = as.numeric(P), POS = as.numeric(POS_hg38),
           BETA = as.numeric(BETA), SE = as.numeric(SE))]

mine <- mine[POS %between% XLIM & is.finite(P)]
his  <- his[POS %between% XLIM & is.finite(P)]

# Attach class for his traits, via the regenie-name key.
tmap_u <- unique(tmap[, .(trait = regenie_trait, class_label, subclass_label,
                          mgp_trait)])
his <- merge(his, tmap_u, by = "trait", all.x = TRUE)
his[trait == "NeuronGliaRatio",
    `:=`(class_label = "NA", subclass_label = "Neuron/glia", mgp_trait = "NA")]
his[, cls := fifelse(grepl("GABA", class_label), "Inhibitory",
             fifelse(grepl("Glutam", class_label), "Excitatory",
             fifelse(grepl("Non-neuronal", class_label), "Glia", "NA")))]

gene_band <- annotate("rect", xmin = GENE[1], xmax = GENE[2],
                      ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.55)
gene_lab <- annotate("text", x = mean(GENE), y = Inf, label = "TMEM106B",
                     vjust = 1.4, size = 3, fontface = "italic",
                     colour = "grey30")

# Panels a/b are SST only. The full trait tables are kept for panel c, which
# still compares every mapped class at his study-wide lead (Lamp5_5), not the
# SST-only peak.
mine_sst <- mine[trait == "SST"]
his_sst  <- his[mgp_trait == "SST"]
stopifnot(nrow(mine_sst) > 0, uniqueN(his_sst$trait) >= 2)

# Shared y-axis so the two SST scans can be read against each other.
ymax <- max(-log10(c(mine_sst$P, his_sst$P)), na.rm = TRUE) * 1.12

# Draw weaker points first so the strongest Sst supertype is not buried.
his_sst <- his_sst[order(-P)]
sst_levels <- his_sst[, .SD[which.min(P)], by = trait][order(P), trait]
his_sst[, trait := factor(trait, levels = sst_levels)]

# A qualitative palette with enough distinct colours for 8 Sst types.
# Sst Chodl is the documented dissenter; give it a colour that stands apart.
SST_COL <- c(
  Sst_25              = "#2E6F8E",
  Sst_3               = "#4C9A8A",
  Sst_11              = "#7FA8BD",
  Sst_19              = "#5B7C99",
  Sst_20              = "#8FBC8F",
  Sst_23              = "#A67C52",
  Sst_12              = "#C4A35A",
  Sst_Chodl_3_SEAAD   = "#C1553B"
)
# Any unexpected Sst name still gets a colour rather than being dropped.
missing_sst <- setdiff(levels(his_sst$trait), names(SST_COL))
if (length(missing_sst)) {
  extra <- grDevices::hcl.colors(length(missing_sst), "Dark 3")
  names(extra) <- missing_sst
  SST_COL <- c(SST_COL, extra)
}

# --- (a) my SST only --------------------------------------------------------
my_sst_best <- mine_sst[which.min(P)]
pa <- ggplot(mine_sst, aes(POS, -log10(P))) +
  gene_band + gene_lab +
  geom_point(colour = PIPE_COL[["mine"]], size = 0.85, alpha = 0.7) +
  geom_hline(yintercept = -log10(GW_LINE), linetype = "22",
             colour = "grey35", linewidth = 0.4) +
  geom_text_repel(data = my_sst_best,
                  aes(label = sprintf("SST  P=%.2g", P)),
                  size = 3.0, fontface = "bold", colour = PIPE_COL[["mine"]],
                  nudge_y = 1.2, min.segment.length = 0) +
  scale_x_continuous(labels = function(x) sprintf("%.2f Mb", x / 1e6)) +
  coord_cartesian(ylim = c(0, ymax)) +
  labs(x = NULL, y = expression(-log[10](P)),
       title = "a  My SST class",
       subtitle = paste("One phenotype, the average of every Sst type.",
                        "Nothing in the window is significant.")) +
  theme_cmp(10)

# --- (b) his Sst-family supertypes only -------------------------------------
his_sst_best <- his_sst[, .SD[which.min(P)], by = trait]
pb <- ggplot(his_sst, aes(POS, -log10(P), colour = trait)) +
  gene_band + gene_lab +
  geom_point(size = 0.85, alpha = 0.75) +
  geom_hline(yintercept = -log10(GW_LINE), linetype = "22",
             colour = "grey35", linewidth = 0.4) +
  geom_hline(yintercept = -log10(MEFF_LINE), linetype = "solid",
             colour = "grey20", linewidth = 0.4) +
  geom_text_repel(data = his_sst_best[P < 5e-8],
                  aes(label = sprintf("%s  P=%.1e", trait, P)),
                  size = 2.6, fontface = "bold", max.overlaps = 12,
                  box.padding = 0.35, show.legend = FALSE) +
  scale_colour_manual(values = SST_COL, name = NULL, drop = FALSE) +
  scale_x_continuous(labels = function(x) sprintf("%.2f Mb", x / 1e6)) +
  coord_cartesian(ylim = c(0, ymax)) +
  labs(x = "chr7 position (hg38)", y = expression(-log[10](P)),
       title = "b  His Sst supertypes (same window, same y-axis)",
       subtitle = paste("Dashed line 5e-8, solid line his corrected",
                        "threshold 2.94e-09. Sst_Chodl is the dissenter.")) +
  theme_cmp(10) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = rel(0.85))) +
  guides(colour = guide_legend(nrow = 2, override.aes = list(size = 2.4,
                                                             alpha = 1)))

# --- (c) what the two phenotypes make of the same allele --------------------
# At his study-wide lead (Lamp5_5), the z of every mapped supertype against the
# z my class-level scan reported at that same variant. Joined on position AND
# alleles so a shared position with a different variant is not treated as a
# match.
lead <- his[which.min(P), variant]
lead_row <- his[variant == lead][1]

his_lead <- his[variant == lead, .(trait, cls, mgp_trait, z = BETA / SE)]
my_lead <- mine[POS == lead_row$POS & toupper(REF) == toupper(lead_row$REF) &
                toupper(ALT) == toupper(lead_row$ALT),
                .(mgp_trait = trait, z_mine = BETA / SE, p_mine = P)]
stopifnot(nrow(my_lead) > 0)
cmp <- merge(his_lead[mgp_trait != "NA"], my_lead, by = "mgp_trait")

# Split by compartment, because that is what the panel is about.
cmp[, compartment := fifelse(cls == "Glia", "glia", "neuronal")]
cmp[, mgp_f := factor(mgp_trait,
                      levels = unique(cmp[order(compartment, -z_mine),
                                          mgp_trait]))]
my_pts <- unique(cmp[, .(mgp_f, compartment, z_mine, p_mine)])

pc <- ggplot(cmp, aes(z, mgp_f)) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.4) +
  geom_segment(data = my_pts, aes(x = z_mine, xend = 0, y = mgp_f,
                                  yend = mgp_f),
               colour = PIPE_COL[["mine"]], alpha = 0.3, linewidth = 1.6) +
  geom_point(aes(colour = cls), size = 2.4, alpha = 0.85,
             position = position_jitter(height = 0.10, width = 0)) +
  geom_point(data = my_pts, aes(x = z_mine, y = mgp_f),
             shape = 23, size = 3.6, stroke = 1.1,
             fill = PIPE_COL[["mine"]], colour = "black") +
  geom_text(data = my_pts[mgp_f == "SST"],
            aes(x = z_mine, y = mgp_f,
                label = sprintf("my SST: z=%.2f, P=%.2f", z_mine, p_mine)),
            vjust = -1.4, hjust = -0.05, size = 2.9, fontface = "bold",
            colour = PIPE_COL[["mine"]]) +
  facet_grid(compartment ~ ., scales = "free_y", space = "free_y") +
  scale_colour_manual(values = CLASS_COL, name = "his supertype, class") +
  labs(x = sprintf("z at his lead variant %s", lead), y = NULL,
       title = "c  Two phenotypes, one allele, two different readings",
       subtitle = paste0(
         "Circles are his supertypes, the red diamond my class-level scan.\n",
         "Glia: my scan reads every glial class strongly negative while his ",
         "supertypes sit at zero.\nThat is whole-composition closure against ",
         "CLR closed within compartment.\n",
         "SST: his 8 Sst supertypes average z +5.95; my SST class reads 0.53.")) +
  theme_cmp(10) +
  theme(legend.position = "bottom",
        panel.grid.major.y = element_line(colour = "grey94"))

p <- ((pa / pb) | pc) + plot_layout(widths = c(1.15, 1)) +
  plot_annotation(
    title = "Stage 2: SST at TMEM106B, class average against Sst supertypes",
    caption = paste(
      "Panels a and b show SST only, on a shared y-axis.",
      "Panel c is every mapped class at his study-wide lead, not the SST peak.",
      "\nMine: results/meta_analysis_15cohorts_hg19_v2/*.annotated.tsv,",
      "effect oriented to ALT, positions carried to hg38 by the offset in",
      "data/coord_offsets.tsv.",
      "\nHis: ctpgwas/results/step2/all_shrunk, pooled arm.",
      "Sst family from celltype-composition/refs/taxonomy_DFC_2026.tsv."),
    theme = theme_cmp(12))

save_fig(p, "fig3_stage2_traits_tmem106b.png", width = 15, height = 10)

# Numbers the summary quotes, recomputed here so the two cannot disagree.
cat("\n  lead variant:", lead, "at",
    format(lead_row$POS, big.mark = ","), "\n")
out <- cmp[, .(n_supertypes = .N, n_positive = sum(z > 0),
               n_negative = sum(z < 0), his_mean_z = round(mean(z), 2),
               my_z = round(z_mine[1], 2), my_p = signif(p_mine[1], 3)),
           by = .(compartment, mgp_trait)][order(compartment, mgp_trait)]
print(out)
fwrite(out, file.path(DATA_DIR, "fig3_lead_comparison.tsv"), sep = "\t")
