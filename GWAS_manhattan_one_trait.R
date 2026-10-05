#!/usr/bin/env Rscript

# ============================================================
# Manhattan plot for ONE GWAS trait
#
# Input:
#   1. MLMA file
#   2. Plot title
#   3. Output directory
#
# Thresholds:
#   Significant: P < 1.10e-8
#   Suggestive : P < 2.20e-7
#
# Plot:
#   - Chromosomes: black / grey
#   - Significant threshold: red line
#   - Suggestive threshold: blue line
#   - Significant SNPs: red
#   - Significant SNP size unchanged
# ============================================================


# ============================================================
# 0. R library
# ============================================================

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


# ============================================================
# 1. Arguments
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {

  stop(
    paste0(
      "Usage:\n",
      "Rscript GWAS_manhattan_one_trait.R ",
      "<MLMA_file> <plot_title> <output_dir>\n"
    )
  )
}

mlma_file <- args[1]
plot_title <- args[2]
output_dir <- args[3]


# ============================================================
# 2. Check input
# ============================================================

if (!file.exists(mlma_file)) {
  stop("MLMA file does not exist: ", mlma_file)
}

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


cat("============================================\n")
cat("MLMA file :", mlma_file, "\n")
cat("Title     :", plot_title, "\n")
cat("Output    :", output_dir, "\n")
cat("============================================\n")


# ============================================================
# 3. Thresholds
# ============================================================

threshold_sig <- 1.10e-8

threshold_suggestive <- 2.20e-7


cat(
  "Significant threshold:",
  format(threshold_sig, scientific = TRUE),
  "\n"
)

cat(
  "Suggestive threshold:",
  format(threshold_suggestive, scientific = TRUE),
  "\n"
)

cat(
  "-log10 significant:",
  -log10(threshold_sig),
  "\n"
)

cat(
  "-log10 suggestive:",
  -log10(threshold_suggestive),
  "\n"
)


# ============================================================
# 4. Read MLMA
# ============================================================

df <- fread(
  mlma_file,
  header = TRUE
)


required_cols <- c(
  "Chr",
  "SNP",
  "bp",
  "p"
)


if (!all(required_cols %in% colnames(df))) {

  stop(
    "Required columns missing from ",
    mlma_file,
    "\nColumns found: ",
    paste(
      colnames(df),
      collapse = ", "
    )
  )
}


# ============================================================
# 5. Remove invalid P values
# ============================================================

df <- df[
  !is.na(p) &
  is.finite(p) &
  p > 0 &
  p <= 1
]


if (nrow(df) == 0) {
  stop("No valid SNPs found.")
}


cat(
  "Valid SNPs:",
  nrow(df),
  "\n"
)

cat(
  "Minimum P:",
  format(
    min(df$p),
    scientific = TRUE
  ),
  "\n"
)


# ============================================================
# 6. Prepare CMplot data
# ============================================================

gwas <- data.frame(

  SNP = df$SNP,

  CHR = suppressWarnings(
    as.integer(df$Chr)
  ),

  BP = df$bp,

  P = df$p
)


gwas <- gwas[
  !is.na(gwas$CHR) &
  !is.na(gwas$BP) &
  !is.na(gwas$P),
]


if (nrow(gwas) == 0) {
  stop("No valid SNPs remain for plotting.")
}


# ============================================================
# 7. Significant SNPs
#
# Only P < 1.10e-8 are highlighted red.
#
# Suggestive SNPs are NOT recoloured.
# ============================================================

sig_snps <- gwas$SNP[
  gwas$P < threshold_sig
]


cat(
  "Significant SNPs (P < 1.10e-8):",
  length(sig_snps),
  "\n"
)


cat(
  "Suggestive SNPs (P < 2.20e-7):",
  sum(gwas$P < threshold_suggestive),
  "\n"
)


# If there are no significant SNPs,
# CMplot should receive NULL

if (length(sig_snps) == 0) {
  sig_snps <- NULL
}


# ============================================================
# 8. Manhattan plot
# ============================================================

old_wd <- getwd()

setwd(output_dir)


CMplot(

  gwas,

  # ----------------------------------------------------------
  # Manhattan plot
  # ----------------------------------------------------------

  plot.type = "m",

  type = "p",

  LOG10 = TRUE,


  # ----------------------------------------------------------
  # Black / grey chromosome colours
  # ----------------------------------------------------------

  col = c(
    "black",
    "grey60"
  ),


  # ----------------------------------------------------------
  # SNP size
  # ----------------------------------------------------------

  cex = 0.4,

  pch = 19,


  # ----------------------------------------------------------
  # Thresholds
  #
  # First:
  #   significant 1.10e-8 -> red
  #
  # Second:
  #   suggestive 2.20e-7 -> blue
  # ----------------------------------------------------------

  threshold = c(
    threshold_sig,
    threshold_suggestive
  ),

  threshold.col = c(
    "red",
    "blue"
  ),

  threshold.lty = c(
    2,
    2
  ),

  threshold.lwd = c(
    1.5,
    1.5
  ),


  # ----------------------------------------------------------
  # Do NOT let CMplot automatically amplify all SNPs
  # crossing the lower threshold.
  #
  # We highlight only genome-wide significant SNPs manually.
  # ----------------------------------------------------------

 # ============================================================
# Significant SNP amplification
# ============================================================

amplify = TRUE,

signal.col = c(
  "red",
   "blue"
),

signal.cex = c(
  0.8,
  0.4
),

signal.pch = c(
  19,
  19
),


  # ----------------------------------------------------------
  # No chromosome density strip
  # ----------------------------------------------------------

  chr.den.col = NULL,


  # ----------------------------------------------------------
  # Title
  # ----------------------------------------------------------

  main = plot_title,

  main.cex = 2.2,

  main.font = 2,


  # ----------------------------------------------------------
  # Output
  # ----------------------------------------------------------

  file = "png",

  file.name = plot_title,

  file.output = TRUE,

  dpi = 300,

  width = 14,

  height = 6,

  verbose = TRUE
)


setwd(old_wd)


# ============================================================
# 9. Finish
# ============================================================

cat("\n============================================\n")
cat("Finished:", plot_title, "\n")
cat("Input   :", mlma_file, "\n")
cat("Output  :", output_dir, "\n")
cat(
  "Top P   :",
  format(
    min(gwas$P),
    scientific = TRUE
  ),
  "\n"
)
cat("============================================\n")