# UpSet plot: which other ROSMAP snRNA-seq datasets also measured the CUIMC2
# donors. Universe = the 237 CUIMC2 persons with >= 1 QC-passing nucleus
# (cuimc2_overlap_data.py, cleaned snRNA vintage 2026-10-13, joined on
# registry person_id). Every column is an EXCLUSIVE membership pattern, so the
# column bars sum to 237; the column with no other dataset (teal, leftmost)
# is the number of CUIMC2 donors not measured by any other ROSMAP snRNA-seq.
#   top:    persons per exclusive pattern
#   left:   persons each dataset shares with CUIMC2 (inclusive)
#   matrix: filled dot = person present in that dataset
# Datasets sharing nobody with CUIMC2 (PsychAD_ROSMAP) are left out of the
# matrix and named in the console output.
# No UpSet package in r_env, so the three panels are drawn with ggplot2 and
# assembled with patchwork.
# Run: /home/zhoux156/miniforge3/envs/r_env/bin/Rscript cuimc2_upset.R
source("/project/rrg-shreejoy/zhoux156/external-Xiaolin/ctp_reliability_paper/common/theme_cns.R")

HERE <- "/project/rrg-shreejoy/zhoux156/external-Xiaolin/ctp_reliability_paper/cohort_overlap"
m <- fread(file.path(HERE, "data/cuimc2_membership.tsv"))
others <- setdiff(names(m), c("person_id", "n_nuclei_CUIMC2"))
sizes <- sapply(others, function(d) sum(m[[d]]))
message("left out (0 shared): ", paste(others[sizes == 0], collapse = ", "))
others <- others[sizes > 0]
sizes <- sizes[others]

# display labels; Multiome regions kept apart because they are different tissue
LAB <- c(CUIMC2 = "CUIMC2", Mathys = "Mathys", Green = "Green",
         Multiome2025_DLPFC = "Diverse multiome DLPFC",
         Multiome2025_STG = "Diverse multiome STG",
         Cain_2023 = "Cain 2023", Zhou_2020 = "Zhou 2020", Mathys_2019 = "Mathys 2019")
LAB["CUIMC2"] <- sprintf("CUIMC2 (n = %d)", nrow(m))
rows <- c("CUIMC2", names(sort(sizes, decreasing = TRUE)))    # top -> bottom
row_lab <- factor(LAB[rows], levels = LAB[rev(rows)])

# exclusive patterns
m[, pattern := apply(.SD, 1, paste, collapse = ""), .SDcols = others]
pat <- m[, .(n = .N), by = pattern]
pat[, n_sets := nchar(gsub("0", "", pattern))]
setorder(pat, -n, n_sets)
pat[, col := factor(seq_len(.N))]
pat[, new := n_sets == 0]
stopifnot(sum(pat$n) == nrow(m))

grid <- pat[, {
  inn <- strsplit(pattern, "")[[1]] == "1"
  .(dataset = c("CUIMC2", others), present = c(TRUE, inn))
}, by = .(col, new)]
grid[, row := factor(LAB[dataset], levels = levels(row_lab))]
seg <- grid[present == TRUE, .(lo = min(as.integer(row)), hi = max(as.integer(row))), by = col]

BAR <- GREYS[3]; NEW <- ACCENT; OFF <- "#E3E3E3"

p_top <- ggplot(pat, aes(col, n, fill = new)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n), vjust = -0.4, size = fontsize / .pt) +
  scale_fill_manual(values = c(`FALSE` = BAR, `TRUE` = NEW), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "CUIMC2 donors\n(exclusive pattern)") +
  theme_cns + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                    axis.line.x = element_blank())

p_mat <- ggplot(grid, aes(col, row)) +
  geom_point(data = grid[present == FALSE], colour = OFF, size = 2.2) +
  geom_segment(data = seg, aes(x = col, xend = col, y = lo, yend = hi),
               inherit.aes = FALSE, linewidth = 0.6, colour = "black") +
  geom_point(data = grid[present == TRUE], aes(colour = new), size = 2.2) +
  scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = NEW), guide = "none") +
  scale_y_discrete(drop = FALSE) +
  labs(x = NULL, y = NULL) +
  theme_cns + theme(axis.text.x = element_blank(), axis.ticks = element_blank(),
                    axis.line = element_blank())

# set sizes: shared with CUIMC2 (the CUIMC2 row itself carries n in its label)
set_dt <- data.table(row = factor(LAB[rows[-1]], levels = levels(row_lab)),
                     n = unname(sizes[rows[-1]]))
p_left <- ggplot(set_dt, aes(n, row)) +
  geom_col(width = 0.6, fill = GREYS[2], orientation = "y") +
  geom_text(aes(label = n), hjust = 1.2, size = fontsize / .pt) +
  scale_x_reverse(expand = expansion(mult = c(0.3, 0))) +
  scale_y_discrete(drop = FALSE) +
  labs(x = "Donors shared\nwith CUIMC2", y = NULL) +
  theme_cns + theme(axis.text.y = element_blank(), axis.line.y = element_blank(),
                    axis.ticks.y = element_blank())

fig <- (plot_spacer() + p_top + p_left + p_mat) +
  plot_layout(ncol = 2, widths = c(1, 3.2), heights = c(1.3, 1))
dir.create(file.path(HERE, "figures"), showWarnings = FALSE)
save_fig(fig, file.path(HERE, "figures/cuimc2_donor_overlap_upset"), height_in = 3.2, width_mm = 150)
print(pat[, .(pattern = sapply(strsplit(pattern, ""), function(z) {
  s <- others[z == "1"]; if (length(s)) paste(s, collapse = " + ") else "CUIMC2 only" }), n)])
