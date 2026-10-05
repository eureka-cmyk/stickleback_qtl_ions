#!/usr/bin/env Rscript

# Generate one Manhattan plot using the established EDA plot settings.
# Usage: Rscript GWAS_compare_one_trait.R <MLMA_file> <plot_title> <output_dir>

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
  stop("Usage: Rscript GWAS_compare_one_trait.R <MLMA_file> <plot_title> <output_dir>")
}

mlma_file <- args[1]
plot_title <- args[2]
output_dir <- args[3]

if (!file.exists(mlma_file)) {
  stop("MLMA file does not exist: ", mlma_file)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

threshold_sig <- 1.10e-8
threshold_suggestive <- 2.20e-7

cat("============================================\n")
cat("MLMA file :", mlma_file, "\n")
cat("Title     :", plot_title, "\n")
cat("Output    :", output_dir, "\n")
cat("============================================\n")

required_cols <- c("Chr", "SNP", "bp", "p")
header <- names(fread(mlma_file, nrows = 0))
if (!all(required_cols %in% header)) {
  stop(
    "Required columns missing from ", mlma_file,
    "\nColumns found: ", paste(header, collapse = ", ")
  )
}

mlma <- fread(mlma_file, select = required_cols)
setnames(mlma, old = c("Chr", "bp", "p"), new = c("CHR", "BP", "P"))
mlma[, CHR := suppressWarnings(as.integer(CHR))]
mlma[, BP := suppressWarnings(as.numeric(BP))]
mlma[, P := suppressWarnings(as.numeric(P))]

gwas <- mlma[
  !is.na(CHR) & is.finite(BP) & BP >= 0 &
    !is.na(P) & P > 0 & P <= 1 & !is.na(SNP),
  .(SNP, CHR, BP, P)
]
setcolorder(gwas, c("SNP", "CHR", "BP", "P"))
rm(mlma)
invisible(gc())
if (nrow(gwas) == 0) {
  stop("No valid SNPs remain for plotting. Check that Chr is numeric.")
}

cat("Valid SNPs:", nrow(gwas), "\n")
cat("Minimum P :", format(min(gwas$P), scientific = TRUE), "\n")
cat("P < 1.10e-8:", sum(gwas$P < threshold_sig), "\n")
cat("P < 2.20e-7:", sum(gwas$P < threshold_suggestive), "\n")

old_wd <- getwd()
setwd(output_dir)

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
  file.name = plot_title,
  file.output = TRUE,
  dpi = 300,
  width = 14,
  height = 6,
  verbose = TRUE
)

setwd(old_wd)
cat("Finished:", plot_title, "\n")
