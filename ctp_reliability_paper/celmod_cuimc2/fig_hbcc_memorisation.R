# Supplement-style diagnostic: NIMH_HBCC bulk x Catching HBCC snRNA CelMod leg (89 paired
# donors, exploratory) memorises its training donors and does not generalise.
#   a  per supertype: in-sample r (model scored on the donors it was fitted to) vs
#      out-of-fold r (5-fold, donors held out) for this pair; points on the dashed
#      identity line would mean no overfitting
#   b  the same two distributions for every pair that has a leg, ordered by number of
#      paired training donors, so the gap can be compared with better-powered pairs
# Inputs (read-only): per-supertype r_oof / r_insample from Shreejoy's
#   results_2026-10-21c/deconv/oof_vs_insample.csv (existing legs) and the two legs built
#   here (scratch). Not a paper figure; for the group discussion of whether the pair is used.
# Run: /home/zhoux156/miniforge3/envs/r_env/bin/Rscript fig_hbcc_memorisation.R
HERE <- "/project/rrg-shreejoy/zhoux156/external-Xiaolin/ctp_reliability_paper"
source(file.path(HERE, "common", "theme_cns.R"))
FIG <- file.path(HERE, "celmod_cuimc2", "figures"); dir.create(FIG, showWarnings = FALSE)
SC <- "/scratch/zhoux156/ctp_reliability_paper"
d <- rbind(
  fread("/scratch/shreejoy/composition/results_2026-10-21c/deconv/oof_vs_insample.csv"),
  fread(file.path(SC, "celmod_cuimc2", "oof_vs_insample__CUIMC2.csv")),
  fread(file.path(SC, "celmod_catching_hbcc", "oof_vs_insample__Catching_2026_HBCC.csv")))
d <- d[level == "supertype"]
# one pair label per bulk x snRNA leg; keep the legs used in Fig 4b plus the two new ones
LEGS <- list(c("ROSMAP","Green"), c("ROSMAP","Mathys"), c("ROSMAP","Multiome2025_DLPFC"),
             c("NIMH_HBCC","PsychAD_HBCC"), c("NABEC","Catching_2026_NABEC"),
             c("ROSMAP","CUIMC2"), c("NIMH_HBCC","Catching_2026_HBCC"))
d <- rbindlist(lapply(LEGS, function(l) d[cohort == l[1] & trainer == l[2]]))
SHORT <- c(Multiome2025_DLPFC = "Diverse_multiome", Catching_2026_NABEC = "Catching", Catching_2026_HBCC = "Catching")
d[, pair := paste0(cohort, " ×\n", gsub("_", "_\n", ifelse(trainer %in% names(SHORT), SHORT[trainer], trainer)))]
d[, pair := ifelse(trainer == "Catching_2026_HBCC", paste0(pair, " HBCC"), pair)]
d[, is_hbcc := trainer == "Catching_2026_HBCC"]
d <- d[is.finite(r_oof) & is.finite(r_insample)]
stopifnot(d[is_hbcc == TRUE, uniqueN(n_train)] == 1, d[is_hbcc == TRUE, n_train[1]] == 89)

smry <- d[, .(n = n_train[1], med_in = median(r_insample), med_oof = median(r_oof),
              n_ct = .N), by = .(pair, is_hbcc)]
smry[, gap := med_in - med_oof]
setorder(smry, -n)
cat("median in-sample vs out-of-fold r (supertype), by pair:\n"); print(smry)
fwrite(smry, file.path(HERE, "celmod_cuimc2", "figures", "hbcc_memorisation_summary.tsv"), sep = "\t")

h <- d[is_hbcc == TRUE]
mh <- smry[is_hbcc == TRUE]
p_a <- ggplot(h, aes(r_insample, r_oof)) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  geom_point(size = 1.1, alpha = 0.7, colour = GREYS[3], stroke = 0) +
  annotate("point", x = mh$med_in, y = mh$med_oof, shape = 23, size = 3, fill = "#E64B35", colour = "black", stroke = 0.4) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.4, size = fontsize / .pt,
           label = "NIMH_HBCC × Catching") +
  coord_cartesian(xlim = c(0.3, 0.8), ylim = c(-0.45, 0.8)) +
  labs(x = "In-sample r (donors the model was fitted on)",
       y = "Out-of-fold r (donors held out)") +
  theme_cns

d[, pair := factor(pair, levels = smry$pair)]
dl <- melt(d, id.vars = c("pair", "is_hbcc", "cell_type", "n_train"),
           measure.vars = c("r_insample", "r_oof"), variable.name = "kind", value.name = "r")
dl[, kind := factor(kind, c("r_insample", "r_oof"), c("In-sample", "Held-out (5-fold)"))]
lab <- smry[, .(pair, lab = sprintf("n = %d", n))]
dl <- merge(dl, lab, by = "pair"); dl[, pair := factor(pair, levels = smry$pair)]
dl[, kg := interaction(kind, is_hbcc)]
p_b <- ggplot(dl, aes(pair, r, fill = kind)) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_boxplot(aes(alpha = is_hbcc), outlier.size = 0.3, outlier.alpha = 0.5, linewidth = 0.3,
               width = 0.7, position = position_dodge(0.78)) +
  scale_fill_manual(values = c("In-sample" = GREYS[1], "Held-out (5-fold)" = ACCENT), name = NULL) +
  scale_alpha_manual(values = c("FALSE" = 0.55, "TRUE" = 1), guide = "none") +
  scale_x_discrete(labels = setNames(paste0(levels(dl$pair), "\nn = ", smry$n), levels(dl$pair))) +
  labs(x = NULL, y = "r, deconvolved vs snRNA-seq (supertypes)") +
  theme_cns + theme(legend.position = "top", legend.margin = margin(0, 0, 0, 0),
                    axis.text.x = element_text(size = fontsize * 0.65, lineheight = 0.9))
fig <- wrap_plots(A = p_a, B = p_b, nrow = 1) + plot_layout(widths = c(1, 1.5)) +
  plot_annotation(tag_levels = "a")
save_fig(fig, file.path(FIG, "hbcc_catching_memorisation"), height_in = 3.6)
