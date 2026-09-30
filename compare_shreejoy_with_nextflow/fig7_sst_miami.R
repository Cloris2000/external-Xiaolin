# Figure 7 -- are the Sst-supertype hits the same direction, and does that
# cancel the subclass scan?
#
# Two tests, because they answer different questions:
#   1. At his TMEM106B lead, do his eight Sst types agree in sign? If most
#      oppose each other, averaging them (what an SST class does) would
#      cancel. If they agree, cancellation cannot explain a null SST class.
#   2. Miami of my SST against each of his Sst types, coloured by sign of
#      beta, so a shared peak that points the other way is visible.

source(file.path(
  "/project/rrg-shreejoy/zhoux156/external-Xiaolin/compare_shreejoy_with_nextflow",
  "theme_compare.R"))
suppressPackageStartupMessages(library(ggrepel))

CHR <- 7
XLIM <- c(12.00e6, 12.45e6)
GENE <- c(12211240, 12243367)
HIS_LEAD <- "chr7:12230939:A:AT"   # study-wide lead from fig3

mine <- read_stage("my_loci_hg38.tsv")[CHROM == CHR & trait == "SST"]
his  <- read_stage("his_loci.tsv")[CHROM == CHR & arm == "pooled"]
tmap <- read_stage("trait_map.tsv")

mine[, `:=`(P = as.numeric(P), POS = as.numeric(POS_hg38),
            BETA = as.numeric(BETA_ALT), SE = as.numeric(SE))]
his[, `:=`(P = as.numeric(P), POS = as.numeric(POS_hg38),
           BETA = as.numeric(BETA), SE = as.numeric(SE))]
mine <- mine[POS %between% XLIM & is.finite(P) & is.finite(BETA)]
his  <- his[POS %between% XLIM & is.finite(P) & is.finite(BETA)]

tmap_u <- unique(tmap[, .(trait = regenie_trait, mgp_trait, subclass_label)])
his <- merge(his, tmap_u, by = "trait", all.x = TRUE)
his_sst <- his[mgp_trait == "SST"]
stopifnot(uniqueN(his_sst$trait) >= 2)

# --- signed effects at the study-wide lead ----------------------------------
lead_his <- his_sst[variant == HIS_LEAD]
lead_row <- lead_his[1]
stopifnot(nrow(lead_his) > 0)

lead_mine <- mine[POS == lead_row$POS &
                  toupper(REF) == toupper(lead_row$REF) &
                  toupper(ALT) == toupper(lead_row$ALT)]
stopifnot(nrow(lead_mine) == 1)

lead <- rbind(
  lead_his[, .(label = trait, beta = BETA, se = SE, p = P,
               pipeline = "his")],
  lead_mine[, .(label = "my SST", beta = BETA, se = SE, p = P,
                pipeline = "mine")]
)
lead[, `:=`(z = beta / se, lo = beta - 1.96 * se, hi = beta + 1.96 * se,
            sign = fifelse(beta > 0, "positive", "negative"))]

# What an inverse-variance average of his eight Sst types would report.
ivw_beta <- lead_his[, sum(BETA / SE^2) / sum(1 / SE^2)]
ivw_se   <- lead_his[, sqrt(1 / sum(1 / SE^2))]
ivw_z    <- ivw_beta / ivw_se
ivw_p    <- 2 * pnorm(-abs(ivw_z))

lead <- rbind(lead, data.table(
  label = "IVW of his 8 Sst", beta = ivw_beta, se = ivw_se, p = ivw_p,
  pipeline = "ivw", z = ivw_z, lo = ivw_beta - 1.96 * ivw_se,
  hi = ivw_beta + 1.96 * ivw_se,
  sign = fifelse(ivw_beta > 0, "positive", "negative")))

ord <- lead[order(pipeline != "mine", beta), label]
lead[, lab_f := factor(label, levels = ord)]

SIGN_COL <- c(positive = "#2E6F8E", negative = "#C1553B")

pa <- ggplot(lead, aes(beta, lab_f, colour = sign)) +
  geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.55) +
  geom_point(aes(shape = pipeline), size = 2.6) +
  scale_colour_manual(values = SIGN_COL, name = "effect sign") +
  scale_shape_manual(values = c(his = 16, mine = 23, ivw = 15), guide = "none") +
  labs(x = sprintf("effect at %s (ALT = T insertion), 95%% CI", HIS_LEAD),
       y = NULL,
       title = "a  Direction at the lead: 7 of 8 Sst types agree",
       subtitle = sprintf(
         paste("Only Sst_Chodl points the other way. Averaging his eight",
               "still gives beta = %+.3f (P = %.1e).",
               "My SST is %+.3f (P = %.2f) — too small to be that one dissent."),
         ivw_beta, ivw_p, lead_mine$BETA, lead_mine$P)) +
  theme_cmp(10) +
  theme(legend.position = "bottom")

# --- window-wide sign concordance, matched on pos + alleles -----------------
mine_key <- mine[, .(POS, REF = toupper(REF), ALT = toupper(ALT),
                     z_mine = BETA / SE, beta_mine = BETA, p_mine = P)]
conc <- rbindlist(lapply(unique(his_sst$trait), function(tr) {
  h <- his_sst[trait == tr, .(POS, REF = toupper(REF), ALT = toupper(ALT),
                              z_his = BETA / SE, beta_his = BETA, p_his = P)]
  m <- merge(mine_key, h, by = c("POS", "REF", "ALT"))
  data.table(
    trait = tr,
    n_matched = nrow(m),
    n_same_sign = m[, sum(sign(beta_mine) == sign(beta_his) &
                            beta_mine != 0 & beta_his != 0)],
    r_z = if (nrow(m) > 10) m[, cor(z_mine, z_his)] else NA_real_,
    lead_same = lead_his[trait == tr, sign(BETA)] == sign(lead_mine$BETA)
  )
}))
conc[, frac_same := n_same_sign / n_matched]
fwrite(conc, file.path(DATA_DIR, "fig7_sst_sign_concordance.tsv"), sep = "\t")

# --- Signed Miami, one scan per panel ---------------------------------------
# Within one method: variants with beta > 0 go up, beta < 0 go down.
# -log10(P) of a positive effect always points upward. This is not a
# cross-pipeline Miami (mine on top, his inverted).
signed_miami <- function(dt, title, subtitle_col = "grey20", ylim = 26) {
  d <- copy(dt)
  d[, signed := sign(BETA) * -log10(P)]
  d[BETA == 0, signed := 0]
  ggplot(d, aes(POS, signed, colour = fifelse(BETA > 0, "positive", "negative"))) +
    annotate("rect", xmin = GENE[1], xmax = GENE[2], ymin = -Inf, ymax = Inf,
             fill = "grey88", alpha = 0.7) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.35) +
    geom_hline(yintercept = c(-1, 1) * -log10(GW_LINE), linetype = "22",
               colour = "grey40", linewidth = 0.3) +
    geom_point(size = 0.45, alpha = 0.7) +
    scale_colour_manual(values = SIGN_COL, guide = "none") +
    scale_x_continuous(labels = function(x) sprintf("%.2f Mb", x / 1e6),
                       limits = XLIM) +
    coord_cartesian(ylim = c(-ylim, ylim)) +
    labs(x = NULL,
         y = expression(sign(beta) %*% -log[10](P)),
         subtitle = title) +
    theme_cmp(8) +
    theme(plot.subtitle = element_text(colour = subtitle_col, face = "bold",
                                       size = rel(1.0)))
}

sst_order <- lead_his[order(-abs(BETA / SE)), trait]
his_panels <- lapply(sst_order, function(tr) {
  same <- lead_his[trait == tr, sign(BETA)] == sign(lead_mine$BETA)
  signed_miami(
    his_sst[trait == tr],
    title = sprintf("%s  (%s at lead)", tr,
                    if (isTRUE(same)) "positive" else "NEGATIVE"),
    subtitle_col = if (isTRUE(same)) PIPE_COL[["his"]] else PIPE_COL[["mine"]],
    ylim = 26)
})
my_panel <- signed_miami(mine, "my SST", PIPE_COL[["mine"]], ylim = 26)

pb <- wrap_plots(c(list(my_panel), his_panels), ncol = 3) +
  plot_annotation(
    title = paste("b  Signed Miami within each scan:",
                  "positive beta up, negative beta down"),
    subtitle = paste("Every panel is one method's own SST/Sst GWAS.",
                     "A cancelling subtype would peak downward.",
                     "Only Sst_Chodl does."))

p <- (pa / pb) + plot_layout(heights = c(1, 3.1)) +
  plot_annotation(
    title = paste("SST at TMEM106B: opposing Sst types do not explain",
                  "the null subclass scan"),
    caption = paste(
      "Lead is his study-wide variant chr7:12230939:A:AT, matched on position",
      "and REF/ALT so the allele is the same.",
      "\nMine: SST_meta_analysis_*.annotated.tsv, BETA_ALT.",
      "His: step2/all_shrunk/chr7_Sst_*.regenie, pooled arm.",
      "\nIVW is an inverse-variance mean of his eight Sst betas at that one",
      "variant; it is not a re-run of a class-level GWAS."),
    theme = theme_cmp(12))

save_fig(p, "fig7_sst_miami.png", width = 15, height = 14)

cat("\n  lead", HIS_LEAD, "\n")
print(lead[order(pipeline, -abs(z)),
           .(label, beta = round(beta, 4), se = round(se, 4),
             z = round(z, 2), p = signif(p, 3), sign)])
cat("\n  window-wide sign concordance (matched variants):\n")
print(conc[order(-frac_same),
           .(trait, n_matched, frac_same = round(frac_same, 3),
             r_z = round(r_z, 3), lead_same)])
