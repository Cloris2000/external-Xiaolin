# Figure 10 supplement (DRAFT): how concordant are bulk-deconvolved and snRNA-seq composition?
#   a  Sst_25, same donors: bulk CelMod prediction vs snRNA CLR (both z-scored), per pair
#   b  same donors, disease effects: snRNA beta vs bulk beta across the 86 Core Supertypes
#   c  same donors, per-donor agreement for all 86 Core Supertypes (Pearson r), Sst_25 marked
#   d  permutation null for the c correlation: case labels permuted (shared by both modalities),
#      so the null keeps measurement agreement and between-type correlation; line = observed
#   e  all cohorts (Fig 10 random-effects meta, mostly different donors): snRNA vs bulk beta
# Inputs: data/fig10_concordance_*.tsv (fig10_concordance_data.py), data/fig10_effects_meta.tsv
# Run: /home/zhoux156/miniforge3/envs/r_env/bin/Rscript fig10_concordance.R
args <- commandArgs(trailingOnly = FALSE)
HERE <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)])))
source(file.path(HERE, "..", "common", "theme_cns.R"))
suppressPackageStartupMessages(library(ggrepel))

TYPE <- "Sst_25"
TAX <- "/project/rrg-shreejoy/zhoux156/shreejoy_pipeline/celltype-composition/refs/taxonomy_DFC_2026.tsv"
CLASS_LAB <- c("Neuronal: Glutamatergic" = "Excitatory", "Neuronal: GABAergic" = "Inhibitory",
               "Non-neuronal and Non-neural" = "Non-neuronal")
CLS_COL <- c(Excitatory = "#D1491C", Inhibitory = "#2F4FA8", "Non-neuronal" = "#2E8B3C")
GRP_COL <- c(Control = "#B7B7B7", Case = "#3D3D3D")
DIS_LAB <- c(SCZ = "SCZ", BD = "BP", AD = "AD")
D <- function(f) fread(file.path(HERE, "data", f))

S <- D("fig10_concordance_summary.tsv")
S <- S[order(match(disease, names(DIS_LAB)), -n)]
S[, short := paste0(DIS_LAB[disease], ": ", cohort_label(bulk_cohort), " ×\n",
                    cohort_label(snrna_cohort), " (n = ", n, ")")]
LEV <- S$short
S[, plab := factor(short, levels = LEV)]
setp <- function(x) {
  x[, plab := factor(S$short[match(x$pair, S$pair)], levels = LEV)]
  stopifnot(!anyNA(x$plab)); x
}
# panel b strips also carry the statistics
bl <- S[, paste0(short, "\nρ ", sprintf("%.2f", rho_beta), " (null ", sprintf("%.2f", null_rho_median),
                 "), p_perm ", sprintf("%.2f", p_perm), "\nsign ", round(100 * sign_agree),
                 "%, SD ratio ", sprintf("%.1f", sd_ratio_bulk_over_snrna))]
S[, blab := factor(bl, levels = bl)]
tax <- fread(TAX)
cls <- setNames(CLASS_LAB[tax$class_label], tax$supertype_label)

# ---- a: Sst_25 per donor ----------------------------------------------------------------
dn <- setp(D("fig10_concordance_sst25_donors.tsv"))
dn[, grp := factor(ifelse(case == 1, "Case", "Control"), levels = c("Control", "Case"))]
pa <- ggplot(dn, aes(z_snrna, z_bulk)) +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.25) +
  geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.25) +
  geom_point(aes(colour = grp), size = 0.6, alpha = 0.7) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "black", linewidth = 0.4) +
  geom_text(data = S, aes(x = -Inf, y = Inf, label = sprintf("r = %.2f", sst25_donor_r)),
            hjust = -0.1, vjust = 1.3, size = 2.3, inherit.aes = FALSE) +
  facet_wrap(~plab, nrow = 1) +
  scale_colour_manual(values = GRP_COL, name = NULL) +
  labs(x = "snRNA-seq Sst_25 (CLR, z)", y = "Bulk Sst_25\n(CelMod, z)") +
  theme_cns + theme(strip.text = element_text(size = 6), legend.position = "right")

# ---- c: per-donor r, all Core types -------------------------------------------------------
rr <- setp(D("fig10_concordance_donor_r.tsv"))
rr[, class := factor(cls[cell_type], levels = CLASS_LAB)]
stopifnot(!anyNA(rr$class))
pb <- ggplot(rr, aes(x = plab, y = r)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.25, linetype = "dashed") +
  geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.3, fill = NA) +
  geom_jitter(aes(colour = class), width = 0.15, height = 0, size = 0.5, alpha = 0.7) +
  geom_point(data = rr[cell_type == TYPE], shape = 23, size = 2, fill = "gold", colour = "black") +
  scale_colour_manual(values = CLS_COL, name = NULL) +
  labs(x = NULL, y = "Per-donor r, bulk vs snRNA\n(86 Core Supertypes)",
       subtitle = "◆ Sst_25") +
  theme_cns + theme(axis.text.x = element_text(size = 5, angle = 30, hjust = 1), legend.position = "none",
                    plot.subtitle = element_text(size = 6.5))

# ---- c: disease betas, same donors ------------------------------------------------------
bt <- setp(D("fig10_concordance_betas.tsv"))
bt[, blab := S$blab[match(plab, S$plab)]]
bt[, class := factor(cls[cell_type], levels = CLASS_LAB)]
lim <- range(c(bt$beta_bulk, bt$beta_snrna))
pc <- ggplot(bt, aes(beta_snrna, beta_bulk)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey60", linewidth = 0.3, linetype = "dashed") +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.25) +
  geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.25) +
  geom_point(aes(colour = class), size = 0.8, alpha = 0.8) +
  geom_point(data = bt[cell_type == TYPE], shape = 23, size = 2, fill = "gold", colour = "black") +
  facet_wrap(~blab, nrow = 1) +
  coord_equal(xlim = lim, ylim = lim) +
  scale_colour_manual(values = CLS_COL, name = NULL) +
  labs(x = "snRNA-seq case β (SD)", y = "Bulk case β (SD)",
       subtitle = "SD ratio = SD(bulk β) / SD(snRNA β) across the 86 types; sign = % same sign") +
  theme_cns + theme(strip.text = element_text(size = 5.5), legend.position = "right")

# ---- d: permutation null --------------------------------------------------------------
nl <- setp(D("fig10_concordance_null.tsv"))
pd <- ggplot(nl, aes(null_rho)) +
  geom_histogram(bins = 40, fill = "grey75", colour = NA) +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.25, linetype = "dashed") +
  geom_vline(data = S, aes(xintercept = rho_beta), colour = "#E64B35", linewidth = 0.6) +
  facet_wrap(~plab, nrow = 2) +
  labs(x = "ρ(bulk β, snRNA β) under permuted case labels", y = "Permutations",
       subtitle = "red = observed") +
  theme_cns + theme(strip.text = element_text(size = 5.5), plot.subtitle = element_text(size = 6.5))

# ---- e: meta level, all cohorts ---------------------------------------------------------
m <- D("fig10_effects_meta.tsv")[model == "primary"]
w <- dcast(m, disease + cell_type ~ modality, value.var = c("beta", "se", "fdr"))
w[, class := factor(cls[cell_type], levels = CLASS_LAB)]
w[, dis := factor(DIS_LAB[disease], levels = DIS_LAB)]
st <- w[, .(rho = cor(beta_bulk, beta_snRNA, method = "spearman"),
            agree = mean(sign(beta_bulk) == sign(beta_snRNA)),
            n_both = sum(fdr_bulk < 0.05 & fdr_snRNA < 0.05),
            agree_both = sum(fdr_bulk < 0.05 & fdr_snRNA < 0.05 & sign(beta_bulk) == sign(beta_snRNA)),
            sdr = sd(beta_bulk) / sd(beta_snRNA)), by = dis]
st[, elab := paste0(dis, "\nρ ", sprintf("%.2f", rho), ", sign ", round(100 * agree),
                    "%, SD ratio ", sprintf("%.1f", sdr),
                    "\nFDR<0.05 in both: ", agree_both, "/", n_both, " same sign")]
w[, elab := factor(st$elab[match(dis, st$dis)], levels = st$elab[order(st$dis)])]
lim2 <- range(c(w$beta_bulk - w$se_bulk, w$beta_bulk + w$se_bulk,
                w$beta_snRNA - w$se_snRNA, w$beta_snRNA + w$se_snRNA))
pe <- ggplot(w, aes(beta_snRNA, beta_bulk)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey60", linewidth = 0.3, linetype = "dashed") +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.25) +
  geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.25) +
  geom_errorbar(aes(ymin = beta_bulk - se_bulk, ymax = beta_bulk + se_bulk, colour = class),
                width = 0, linewidth = 0.2, alpha = 0.5) +
  geom_errorbarh(aes(xmin = beta_snRNA - se_snRNA, xmax = beta_snRNA + se_snRNA, colour = class),
                 height = 0, linewidth = 0.2, alpha = 0.5) +
  geom_point(aes(colour = class), size = 0.8) +
  geom_point(data = w[cell_type == TYPE], shape = 23, size = 2, fill = "gold", colour = "black") +
  facet_wrap(~elab, nrow = 1) +
  coord_equal(xlim = lim2, ylim = lim2) +
  scale_colour_manual(values = CLS_COL, name = NULL) +
  labs(x = "snRNA-seq case β, RE meta (SD)", y = "Bulk case β,\nRE meta (SD)",
       subtitle = "All cohorts (Fig 10)") +
  theme_cns + theme(legend.position = "right", plot.subtitle = element_text(size = 6.5),
                    strip.text = element_text(size = 6))

fig <- wrap_elements(full = pa) / wrap_elements(full = pc) /
  (wrap_elements(full = pb) | wrap_elements(full = pd)) / wrap_elements(full = pe) +
  plot_layout(heights = c(1, 1.15, 1.15, 1.15)) + plot_annotation(tag_levels = "a")
# panel tags follow reading order: a Sst_25 donors, b betas (same donors), c donor r, d null, e meta
save_fig(fig, file.path(HERE, "figures", "figure10_supp_concordance_draft"), height_in = 10.5)
