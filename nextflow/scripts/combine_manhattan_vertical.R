#!/usr/bin/env Rscript
# Combine two Manhattan plots vertically
# Usage: Rscript combine_manhattan_vertical.R <plot1.png> <plot2.png> <output.png> <label1> <label2>

suppressPackageStartupMessages({
  library(magick)
})

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript combine_manhattan_vertical.R <plot1.png> <plot2.png> <output.png> [label1] [label2]")
}

plot1_file <- args[1]
plot2_file <- args[2]
output_file <- args[3]
label1 <- if (length(args) >= 4) args[4] else "A"
label2 <- if (length(args) >= 5) args[5] else "B"

cat("============================================================\n")
cat("Combining Manhattan Plots Vertically\n")
cat("============================================================\n")
cat("Plot 1:", plot1_file, "\n")
cat("Plot 2:", plot2_file, "\n")
cat("Output:", output_file, "\n")
cat("Labels:", label1, "|", label2, "\n\n")

# Read images
cat("Reading images...\n")
img1 <- image_read(plot1_file)
img2 <- image_read(plot2_file)

# Get dimensions
info1 <- image_info(img1)
info2 <- image_info(img2)

cat("Plot 1 dimensions:", info1$width, "x", info1$height, "\n")
cat("Plot 2 dimensions:", info2$width, "x", info2$height, "\n")

# Make sure both images have the same width (use the larger width)
target_width <- max(info1$width, info2$width)

if (info1$width != target_width) {
  cat("Resizing plot 1 to width", target_width, "\n")
  img1 <- image_resize(img1, paste0(target_width, "x"))
}

if (info2$width != target_width) {
  cat("Resizing plot 2 to width", target_width, "\n")
  img2 <- image_resize(img2, paste0(target_width, "x"))
}

# Add panel labels (A, B) to each plot
cat("Adding panel labels...\n")
img1 <- image_annotate(img1, label1, 
                       size = 80, 
                       color = "black", 
                       weight = 700,
                       location = "+50+50",
                       font = "Arial")

img2 <- image_annotate(img2, label2, 
                       size = 80, 
                       color = "black", 
                       weight = 700,
                       location = "+50+50",
                       font = "Arial")

# Stack vertically
cat("Stacking plots vertically...\n")
combined <- image_append(c(img1, img2), stack = TRUE)

# Save combined image
cat("Saving combined plot...\n")
image_write(combined, output_file, format = "png")

# Get final dimensions
final_info <- image_info(combined)
cat("\n============================================================\n")
cat("Combined plot saved:", output_file, "\n")
cat("Final dimensions:", final_info$width, "x", final_info$height, "\n")
cat("File size:", file.size(output_file) / 1024^2, "MB\n")
cat("============================================================\n")
