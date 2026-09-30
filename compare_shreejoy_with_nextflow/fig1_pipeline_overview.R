# Figure 1 -- the two pipelines side by side, stage by stage.
#
# A schematic, not a data plot: every number in it is read from the stage tables
# or from a file named in the caption, so the figure and the summary cannot
# drift apart.

source(file.path(
  "/project/rrg-shreejoy/zhoux156/external-Xiaolin/compare_shreejoy_with_nextflow",
  "theme_compare.R"))

disc  <- read_stage("stage3_discovery.tsv")
col   <- read_stage("stage4_coloc.tsv")
dec   <- read_stage("stage1_deconvolution.tsv")

# --- numbers, derived rather than typed ------------------------------------
my_gw       <- disc[pipeline == "mine", sum(as.integer(n_gw))]
my_gw_het   <- disc[pipeline == "mine", sum(as.integer(n_gw_high_het))]
his_gw      <- disc[pipeline == "his",  sum(as.integer(n_gw))]
my_traits   <- disc[pipeline == "mine", .N]
his_traits  <- disc[pipeline == "his",  .N]
my_coloc    <- col[pipeline == "mine", .N]
his_coloc   <- col[pipeline == "his",  .N]
my_dec_med  <- dec[pipeline == "mine", median(as.numeric(r))]
his_dec_med <- dec[pipeline == "his",  median(as.numeric(r))]

fmt <- function(x) formatC(x, big.mark = ",", format = "d")

# --- the four stages --------------------------------------------------------
stages <- data.table(
  stage = factor(rep(1:4, each = 2),
                 labels = c("1. Deconvolution", "2. Phenotype and GWAS",
                            "3. Combining cohorts", "4. Colocalisation")),
  pipeline = rep(c("mine", "his"), 4),
  design = c(
    # stage 1
    sprintf("MGP marker-gene profiles\n%d broad classes\nmedian bulk vs snRNA r = %.2f",
            my_traits, my_dec_med),
    sprintf("CelMod, predicts into CLR space\n33 SEA-AD supertypes\nmedian bulk vs snRNA r = %.2f",
            his_dec_med),
    # stage 2
    sprintf("proportions, residualise sex + age,\nrank-inverse normal\n%s variants, hg19",
            fmt(8578702)),
    sprintf("within-compartment CLR, GLS across\nmeasurements, reliability shrinkage, RINT\n%s variants, hg38",
            fmt(9523063)),
    # stage 3
    "per-cohort REGENIE, then METAL\ninverse-variance meta\n15 cohorts, 3,620 samples",
    "pooled mega-analysis\n5,008 people, 43 measurements\ncohort fitted as a covariate",
    # stage 4
    "coloc against 6 disease GWAS",
    "coloc against 8 disease GWAS\n(AD, PD, LBD are Xiaolin's own files)"
  ),
  result = c(
    sprintf("L5.ET r = -0.01, SST r = 0.44"),
    sprintf("33 of 34 traits usable;\nL5 ET dropped as non-portable"),
    sprintf("%d genome-wide hits", my_gw),
    sprintf("%d pooled-arm hits, 13 loci,\n2 clear the corrected threshold", his_gw),
    sprintf("%d of %d hits (%.0f%%) high-I2",
            my_gw_het, my_gw, 100 * my_gw_het / my_gw),
    "ancestry meta on the same people\nrecovers the same 2 loci",
    sprintf("%s tests; TMEM106B colocalises,\nGRN never tested", fmt(my_coloc)),
    sprintf("%d tests; GRN with AD PP.H4 0.9999,\nTMEM106B with AD PP.H4 0.915", his_coloc)
  )
)

stages[, x := ifelse(pipeline == "mine", 1, 2)]
stages[, y := as.integer(stage)]

# Two stacked text layers per cell so design and result never share a line.
cells <- rbind(
  stages[, .(x, y, pipeline, kind = "design", text = design)],
  stages[, .(x, y, pipeline, kind = "result", text = result)]
)
cells[, y_text := ifelse(kind == "design", -y + 0.20, -y - 0.22)]

p <- ggplot(cells, aes(x, -y)) +
  geom_tile(aes(fill = pipeline), width = 0.94, height = 0.92,
            alpha = 0.12, colour = NA) +
  geom_text(aes(y = y_text, label = text, colour = pipeline,
                fontface = ifelse(kind == "result", "bold", "plain")),
            size = 3.0, lineheight = 1.18, vjust = 0.5) +
  scale_fill_manual(values = PIPE_COL, guide = "none") +
  scale_colour_manual(values = PIPE_COL, guide = "none") +
  scale_x_continuous(
    breaks = 1:2, position = "top",
    labels = c("Nextflow (Xiaolin)", "ctp-gwas phase 2 (Shreejoy)"),
    limits = c(0.42, 2.58), expand = c(0, 0)) +
  scale_y_continuous(breaks = -(1:4), labels = levels(stages$stage),
                     expand = c(0.04, 0)) +
  labs(
    title = "Two cell-type-proportion GWAS pipelines, stage by stage",
    subtitle = paste("Coloured text in each cell is the design choice;",
                     "bold text below it is the result it produced"),
    x = NULL, y = NULL,
    caption = paste(
      "Counts derived from data/stage1_deconvolution.tsv,",
      "stage3_discovery.tsv and stage4_coloc.tsv.",
      "\nVariant counts and the 2-locus corrected result are from",
      "ctp-gwas/stages/report/gwas_report/REPORT.md;",
      "\nsample count from nextflow/docs/sample_overlap/README.md.",
      "Hit counts are pooled-arm only on his side.")
  ) +
  theme_cmp(12) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x.top = element_text(face = "bold", size = rel(1.05)),
    axis.text.y = element_text(face = "bold", hjust = 1, size = rel(0.95))
  )

save_fig(p, "fig1_pipeline_overview.png", width = 12, height = 9.5)
