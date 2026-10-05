#!/usr/bin/env Rscript

# Generate Manhattan and QQ plots for one GCTA MLMA trait.
# Usage: Rscript GWAS_H_W_MAH_one_trait_eda.R <MLMA_file> <plot_title> <output_dir>

.libPaths(c(
  "/gpfs01/home/shxyt1/R/x86_64-pc-linux-gnu-library/4.4",
  .libPaths()
))

options(bitmapType = "cairo")
if (!capabilities("cairo")) {
  stop("This R installation does not support Cairo graphics.")
}

suppressPackageStartupMessages({
  library(data.table)
  library(CMplot)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(paste0(
    "Usage:\n",
    "Rscript GWAS_H_W_MAH_one_trait_eda.R <MLMA_file> <plot_title> <output_dir>\n"
  ))
}

mlma_file <- args[1]
plot_title <- args[2]
output_dir <- args[3]

if (!file.exists(mlma_file)) {
  stop("MLMA file does not exist: ", mlma_file)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("============================================\n")
cat("MLMA file :", mlma_file, "\n")
cat("Title     :", plot_title, "\n")
cat("Output    :", output_dir, "\n")
cat("============================================\n")

threshold_sig <- 1.10e-8
threshold_suggestive <- 2.20e-7

cat("Significant threshold:", format(threshold_sig, scientific = TRUE), "\n")
cat("Suggestive threshold :", format(threshold_suggestive, scientific = TRUE), "\n")

# Read MLMA results and retain valid P values.
df <- fread(mlma_file, header = TRUE)
required_cols <- c("Chr", "SNP", "bp", "p")
if (!all(required_cols %in% colnames(df))) {
  stop(
    "Required columns missing from ", mlma_file,
    "\nColumns found: ", paste(colnames(df), collapse = ", ")
  )
}

df <- df[
  !is.na(p) & is.finite(p) & p > 0 & p <= 1 &
    !is.na(Chr) & !is.na(bp) & !is.na(SNP)
]
if (nrow(df) == 0) {
  stop("No valid SNPs found.")
}

# CMplot requires numeric chromosome IDs.
gwas <- data.frame(
  SNP = df$SNP,
  CHR = suppressWarnings(as.integer(df$Chr)),
  BP = suppressWarnings(as.numeric(df$bp)),
  P = suppressWarnings(as.numeric(df$p))
)
gwas <- gwas[
  !is.na(gwas$CHR) & is.finite(gwas$BP) & !is.na(gwas$P) &
    gwas$BP >= 0 & gwas$P > 0 & gwas$P <= 1,
]
if (nrow(gwas) == 0) {
  stop("No valid SNPs remain for plotting. Check that Chr is numeric.")
}

cat("Valid SNPs:", nrow(gwas), "\n")
cat("Minimum P :", format(min(gwas$P), scientific = TRUE), "\n")
cat("Significant SNPs (P < 1.10e-8):", sum(gwas$P < threshold_sig), "\n")
cat("Suggestive SNPs (P < 2.20e-7):", sum(gwas$P < threshold_suggestive), "\n")

old_wd <- getwd()
setwd(output_dir)

# Manhattan plot: keep the original chromosome colours, thresholds,
# threshold line styles, point sizes, and signal highlighting.
CMplot(
  gwas,
  plot.type = "m",
  type = "p",
  LOG10 = TRUE,
  col = c("black", "grey60"),
  cex = 0.4,
  pch = 19,
  threshold = c(threshold_sig, threshold_suggestive),
  threshold.col = c("red", "blue"),
  threshold.lty = c(2, 2),
  threshold.lwd = c(1.5, 1.5),
  amplify = TRUE,
  signal.col = c("red", "blue"),
  signal.cex = c(0.8, 0.4),
  signal.pch = c(19, 19),
  chr.den.col = NULL,
  main = plot_title,
  main.cex = 2.2,
  main.font = 2,
  file = "png",
  file.name = paste0(plot_title, "_Manhattan"),
  file.output = TRUE,
  dpi = 300,
  width = 14,
  height = 6,
  verbose = TRUE
)

# QQ plot for the same trait.
CMplot(
  gwas,
  plot.type = "q",
  LOG10 = TRUE,
  main = paste0(plot_title, " QQ plot"),
  main.cex = 2.2,
  main.font = 2,
  file = "png",
  file.name = paste0(plot_title, "_QQ"),
  file.output = TRUE,
  dpi = 300,
  width = 6,
  height = 6,
  verbose = TRUE
)

setwd(old_wd)

cat("\n============================================\n")
cat("Finished:", plot_title, "\n")
cat("Input   :", mlma_file, "\n")
cat("Output  :", output_dir, "\n")
cat("Top P   :", format(min(gwas$P), scientific = TRUE), "\n")
cat("============================================\n")
