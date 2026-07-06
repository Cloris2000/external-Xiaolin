#!/usr/bin/env Rscript
# Forest plot for a single SNP across all per-cohort GWAS files.
# Default: top Oligodendrocyte SNP (chr13:40969426) from the 13-cohort
# design-matrix meta-analysis.
#
# Usage:
#   Rscript forest_plot_top_snp.R
#   Rscript forest_plot_top_snp.R --cell_type Microglia --marker chr6:164862615:T:C

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(optparse)
})

option_list <- list(
  make_option("--cell_type", type = "character", default = "Oligodendrocyte"),
  make_option("--marker",    type = "character", default = "chr13:40969426:T:C",
    help = "MarkerName in CHROM:POS:REF:ALT format"),
  make_option("--results_dir", type = "character",
    default = "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results"),
  make_option("--meta_tbl", type = "character", default = "",
    help = "Path to meta-analysis .tbl file (auto-detected if empty)"),
  make_option("--outdir", type = "character",
    default = "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/sn_bulk_meta_similarity_design_matrix/top_hits"),
  make_option("--meta_dir", type = "character",
    default = "/external/rprshnas01/netdata_kcni/stlab/Xiaolin/nextflow/results/meta_analysis_15cohorts"),
  make_option("--sn_tbl", type = "character", default = "",
    help = "Optional path to snRNA-seq meta-analysis .tbl file; if provided, adds an sn-meta row"),
  make_option("--meta_only", action = "store_true", default = FALSE,
    help = "If set, only plot meta-analysis rows (bulk + sn), no per-cohort rows")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

COHORTS <- c("AMP_AD_Mayo", "AMP_AD_Rush",
             "CMC_MSSM", "CMC_PENN", "CMC_PITT", "GTEx", "GVEX",
             "MSBB", "Mayo", "NABEC",
             "NIMH_HBCC_1M", "NIMH_HBCC_Omni5M", "NIMH_HBCC_h650",
             "ROSMAP", "ROSMAP_array")

marker   <- opt$marker
ct       <- opt$cell_type

# parse marker into components
parts  <- strsplit(marker, ":", fixed = TRUE)[[1]]
chrom_num <- as.integer(sub("^chr", "", parts[1]))
pos_bp    <- as.integer(parts[2])
ref_al    <- toupper(parts[3])
alt_al    <- toupper(parts[4])

message("Looking up: ", marker, "  (", ct, ")")

# ---- complement helper ------------------------------------------------------
comp_base <- function(b) c(A="T",T="A",C="G",G="C")[b]

# ---- read one cohort raw_p file, extract SNP --------------------------------
read_cohort_snp <- function(cohort) {
  path <- file.path(opt$results_dir, cohort, "regenie_step2",
                    paste0(cohort, "_", ct, "_step2.regenie.raw_p"))
  if (!file.exists(path)) {
    message("  MISSING: ", cohort)
    return(NULL)
  }
  dt <- fread(path, sep = " ", data.table = TRUE, verbose = FALSE)
  # regenie columns: CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ N TEST BETA SE CHISQ LOG10P EXTRA P
  hit <- dt[CHROM == chrom_num & GENPOS == pos_bp]
  if (nrow(hit) == 0) {
    message("  SNP not found in ", cohort)
    return(data.table(cohort = cohort, beta = NA_real_, se = NA_real_,
                      p = NA_real_, n = NA_real_, found = FALSE,
                      flipped = FALSE, a0 = NA_character_, a1 = NA_character_))
  }
  h <- hit[1]
  a0 <- toupper(h$ALLELE0)   # reference / non-effect allele
  a1 <- toupper(h$ALLELE1)   # effect allele

  # harmonise: align effect allele to meta alt_al
  flip <- FALSE
  if (a1 == alt_al && a0 == ref_al) {
    flip <- FALSE
  } else if (a1 == ref_al && a0 == alt_al) {
    flip <- TRUE
  } else if (a1 == comp_base(alt_al) && a0 == comp_base(ref_al)) {
    flip <- FALSE
  } else if (a1 == comp_base(ref_al) && a0 == comp_base(alt_al)) {
    flip <- TRUE
  } else {
    message("  Allele mismatch in ", cohort, " (a0=", a0, " a1=", a1, ")")
    return(data.table(cohort = cohort, beta = NA_real_, se = NA_real_,
                      p = NA_real_, n = NA_real_, found = FALSE,
                      flipped = FALSE, a0 = a0, a1 = a1))
  }

  beta <- as.numeric(h$BETA) * ifelse(flip, -1, 1)
  se   <- as.numeric(h$SE)
  p    <- as.numeric(h$P)
  n    <- as.numeric(h$N)

  data.table(cohort = cohort, beta = beta, se = se, p = p, n = n,
             found = TRUE, flipped = flip, a0 = a0, a1 = a1)
}

# ---- collect all cohorts ----------------------------------------------------
if (opt$meta_only) {
  cohort_data <- data.table(cohort = factor(character(0), levels = rev(COHORTS)),
                            beta = numeric(0), se = numeric(0), p = numeric(0),
                            n = numeric(0), found = logical(0), flipped = logical(0),
                            a0 = character(0), a1 = character(0))
} else {
  cohort_data <- rbindlist(lapply(COHORTS, read_cohort_snp), fill = TRUE)
  cohort_data[, cohort := factor(cohort, levels = rev(COHORTS))]
}

# ---- read meta-analysis row -------------------------------------------------
# auto-detect meta tbl
if (opt$meta_tbl == "") {
  meta_files <- list.files(opt$meta_dir,
    pattern = paste0("^", ct, "_meta_analysis_.*\\.tbl$"),
    full.names = TRUE)
  meta_files <- meta_files[!grepl("1\\.tbl$", basename(meta_files))]
  opt$meta_tbl <- if (length(meta_files) > 0) meta_files[1] else ""
}

meta_row <- NULL
if (opt$meta_tbl != "" && file.exists(opt$meta_tbl)) {
  meta_dt <- fread(opt$meta_tbl, sep = "\t", data.table = TRUE, verbose = FALSE)
  setnames(meta_dt, tolower(gsub("[^A-Za-z0-9]", "_", names(meta_dt))))
  hit <- meta_dt[markername == marker]
  if (nrow(hit) > 0) {
    h <- hit[1]
    # compute total N from cohorts that have data
    total_n <- sum(cohort_data$n, na.rm = TRUE)
    meta_row <- data.table(
      cohort  = factor("Meta-analysis", levels = c("Meta-analysis", rev(COHORTS))),
      beta    = as.numeric(h$effect),
      se      = as.numeric(h$stderr),
      p       = as.numeric(h$p_value),
      n       = total_n,
      found   = TRUE,
      flipped = FALSE,
      a0      = toupper(h$allele2),
      a1      = toupper(h$allele1)
    )
  }
}

# ---- read snRNA-seq meta row (optional) -------------------------------------
sn_row <- NULL
if (opt$sn_tbl != "" && file.exists(opt$sn_tbl)) {
  sn_dt <- fread(opt$sn_tbl, sep = "\t", data.table = TRUE, verbose = FALSE)
  setnames(sn_dt, tolower(gsub("[^A-Za-z0-9]", "_", names(sn_dt))))
  # lookup by position (MarkerName may differ in format)
  sn_dt[, chrom_num_col := as.integer(sub("^chr","", sapply(strsplit(markername,":"), `[`, 1)))]
  sn_dt[, pos_col       := as.integer(sapply(strsplit(markername,":"), `[`, 2))]
  sn_hit <- sn_dt[chrom_num_col == chrom_num & pos_col == pos_bp]
  if (nrow(sn_hit) > 0) {
    h <- sn_hit[1]
    sn_a1 <- toupper(sapply(strsplit(h$markername, ":"), `[`, 3))  # REF in MarkerName
    sn_a2 <- toupper(sapply(strsplit(h$markername, ":"), `[`, 4))  # ALT in MarkerName
    # harmonise to bulk alt_al
    sn_beta <- as.numeric(h$effect)
    if (!is.na(sn_beta)) {
      if (sn_a2 == alt_al || sn_a2 == comp_base(alt_al)) {
        # already aligned
      } else if (sn_a1 == alt_al || sn_a1 == comp_base(alt_al)) {
        sn_beta <- -sn_beta
      }
    }
    sn_row <- data.table(
      cohort  = factor("snRNA-seq Meta", levels = c("snRNA-seq Meta", "Meta-analysis", rev(COHORTS))),
      beta    = sn_beta,
      se      = as.numeric(h$stderr),
      p       = as.numeric(h$p_value),
      n       = NA_real_,
      found   = TRUE,
      flipped = FALSE,
      a0      = sn_a1,
      a1      = sn_a2
    )
    message("snRNA-seq meta: beta=", round(sn_beta,4), " se=", h$stderr, " p=", h$p_value)
  } else {
    message("SNP not found in sn meta tbl")
  }
}

# ---- combine for plotting ---------------------------------------------------
all_levels <- c(if (!is.null(sn_row)) "snRNA-seq Meta", "Meta-analysis", rev(COHORTS))
cohort_data[, cohort := factor(as.character(cohort), levels = all_levels)]
if (!is.null(meta_row)) {
  meta_row[, cohort := factor(as.character(cohort), levels = all_levels)]
}
if (!is.null(sn_row)) {
  sn_row[, cohort := factor(as.character(cohort), levels = all_levels)]
}
plot_dt <- rbind(
  if (!is.null(sn_row)) sn_row,
  if (!is.null(meta_row)) meta_row,
  cohort_data,
  fill = TRUE
)

plot_dt[, ci_lo   := beta - 1.96 * se]
plot_dt[, ci_hi   := beta + 1.96 * se]
plot_dt[, p_label := ifelse(!is.na(p),
  ifelse(p < 0.001, formatC(p, format = "e", digits = 2),
         sprintf("%.3f", p)),
  "NA")]
plot_dt[, n_label  := ifelse(!is.na(n),
  formatC(as.integer(n), big.mark = ",", format = "d"), "")]
plot_dt[, is_meta  := cohort == "Meta-analysis"]
plot_dt[, is_sn    := cohort == "snRNA-seq Meta"]
plot_dt[, pt_colour := ifelse(is_meta, "#1a7abf", ifelse(is_sn, "#e07b39", "#444444"))]

# ---- plot -------------------------------------------------------------------
xlim_abs <- max(abs(c(plot_dt$ci_lo, plot_dt$ci_hi)), na.rm = TRUE) * 1.15
# annotation columns placed beyond the CI panel via clip = "off"
x_p <- xlim_abs * 1.12   # p-value column x
x_n <- xlim_abs * 1.60   # N column x

# header row at topmost level
top_cohort <- levels(plot_dt$cohort)[nlevels(plot_dt$cohort)]
header_dt  <- data.table(
  cohort  = factor(top_cohort, levels = levels(plot_dt$cohort)),
  x_p     = x_p,
  x_n     = x_n
)

# separator positions (only meaningful when cohort rows are shown)
n_cohort_rows <- if (opt$meta_only) 0L else length(COHORTS)
sep_meta <- if (!is.null(meta_row) && !opt$meta_only) n_cohort_rows + 0.5 else NULL
sep_sn   <- if (!is.null(sn_row) && !is.null(meta_row) && !opt$meta_only) n_cohort_rows + 1.5 else
            if (!is.null(sn_row) && !opt$meta_only) n_cohort_rows + 0.5 else NULL

p <- ggplot(plot_dt, aes(x = beta, y = cohort)) +
  # shaded meta CI band
  {if (!is.null(meta_row))
    annotate("rect",
      xmin = meta_row$beta - 1.96 * meta_row$se,
      xmax = meta_row$beta + 1.96 * meta_row$se,
      ymin = -Inf, ymax = Inf,
      fill = "#1a7abf", alpha = 0.08)
  } +
  geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.4) +
  # per-cohort CI + points
  geom_errorbarh(
    data    = plot_dt[found == TRUE & is_meta == FALSE & is_sn == FALSE],
    aes(xmin = ci_lo, xmax = ci_hi),
    height  = 0.22, colour = "#555555", linewidth = 0.55
  ) +
  geom_point(
    data    = plot_dt[found == TRUE & is_meta == FALSE & is_sn == FALSE],
    colour  = "#555555", size = 2.4
  ) +
  # meta CI + diamond
  {if (!is.null(meta_row))
    geom_errorbarh(
      data    = plot_dt[is_meta == TRUE],
      aes(xmin = ci_lo, xmax = ci_hi),
      height  = 0.32, colour = "#1a7abf", linewidth = 1.1
    )
  } +
  {if (!is.null(meta_row))
    geom_point(data = plot_dt[is_meta == TRUE],
               shape = 18, size = 5.5, colour = "#1a7abf")
  } +
  # snRNA-seq meta CI + diamond
  {if (!is.null(sn_row))
    geom_errorbarh(
      data    = plot_dt[is_sn == TRUE],
      aes(xmin = ci_lo, xmax = ci_hi),
      height  = 0.32, colour = "#e07b39", linewidth = 1.1
    )
  } +
  {if (!is.null(sn_row))
    geom_point(data = plot_dt[is_sn == TRUE],
               shape = 18, size = 5.5, colour = "#e07b39")
  } +
  # dashed separators
  {if (!is.null(sep_meta))
    geom_hline(yintercept = sep_meta, linetype = "dashed", colour = "grey70", linewidth = 0.4)
  } +
  {if (!is.null(sep_sn))
    geom_hline(yintercept = sep_sn, linetype = "dashed", colour = "#e07b39", linewidth = 0.4)
  } +
  # p-value labels (right of panel, clip = off)
  geom_text(
    data = plot_dt[found == TRUE & is_meta == FALSE & is_sn == FALSE],
    aes(x = x_p, label = p_label),
    hjust = 0.5, size = 3.0, colour = "#333333", show.legend = FALSE
  ) +
  {if (!is.null(meta_row))
    geom_text(data = plot_dt[is_meta == TRUE],
              aes(x = x_p, label = p_label),
              hjust = 0.5, size = 3.2, colour = "#1a7abf",
              fontface = "bold", show.legend = FALSE)
  } +
  {if (!is.null(sn_row))
    geom_text(data = plot_dt[is_sn == TRUE],
              aes(x = x_p, label = p_label),
              hjust = 0.5, size = 3.2, colour = "#e07b39",
              fontface = "bold", show.legend = FALSE)
  } +
  # N labels
  geom_text(
    data = plot_dt[found == TRUE & is_meta == FALSE & is_sn == FALSE],
    aes(x = x_n, label = n_label),
    hjust = 0.5, size = 2.8, colour = "grey45", show.legend = FALSE
  ) +
  {if (!is.null(meta_row))
    geom_text(data = plot_dt[is_meta == TRUE],
              aes(x = x_n, label = n_label),
              hjust = 0.5, size = 3.0, colour = "#1a7abf",
              fontface = "bold", show.legend = FALSE)
  } +
  # column headers above topmost row
  geom_text(data = header_dt,
            aes(x = x_p, y = cohort, label = "p-value"),
            hjust = 0.5, vjust = -1.8, size = 3.2,
            fontface = "bold", colour = "grey20", inherit.aes = FALSE) +
  geom_text(data = header_dt,
            aes(x = x_n, y = cohort, label = "N"),
            hjust = 0.5, vjust = -1.8, size = 3.2,
            fontface = "bold", colour = "grey20", inherit.aes = FALSE) +
  scale_x_continuous(
    limits = c(-xlim_abs * 1.02, xlim_abs * 1.02),
    expand = expansion(0),
    name   = paste0("Beta (95% CI)  [effect allele = ", alt_al, "]")
  ) +
  coord_cartesian(clip = "off") +
  labs(
    title    = paste0("Forest plot: ", marker),
    subtitle = paste0("Cell type: ", ct,
                      "  |  Effect allele: ", alt_al,
                      "  |  Reference allele: ", ref_al,
                      "\nEffect = change in ", ct, " proportion per copy of ", alt_al),
    y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.y = element_line(colour = "grey93", linewidth = 0.3),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(colour = "grey88", linewidth = 0.3),
    axis.text.y        = element_text(size = 11),
    plot.title         = element_text(face = "bold", size = 13),
    plot.subtitle      = element_text(size = 10, colour = "grey30"),
    plot.margin        = margin(t = 25, r = 140, b = 10, l = 10, unit = "pt")
  )

safe_marker <- gsub(":", "_", marker)
out_png <- file.path(opt$outdir,
  paste0("forest_", ct, "_", safe_marker, ".png"))
ggsave(out_png, p, width = 10, height = 7, dpi = 300)
message("Saved: ", out_png)

# ---- also save per-cohort table ---------------------------------------------
out_tsv <- file.path(opt$outdir,
  paste0("forest_", ct, "_", safe_marker, "_cohort_data.tsv"))
fwrite(plot_dt[, .(cohort, beta, se, ci_lo, ci_hi, p, n, found, flipped)],
       out_tsv, sep = "\t")
message("Saved table: ", out_tsv)
