#!/usr/bin/env Rscript

# ============================================================
# GWAS plotting for ONE trait
#
# Outputs:
#   1. Manhattan plot: P < 5e-8
#   2. Manhattan plot: Bonferroni threshold
#   3. QQ plot
#   4. Significant SNP table
#   5. Trait-level summary
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

if (length(args) != 2) {
  stop(
    "Usage: Rscript GWAS_H_W_MAH_one_trait.R ",
    "<MLMA_file> <base_dir>"
  )
}

mlma_file <- args[1]
base_dir  <- args[2]

if (!file.exists(mlma_file)) {
  stop("MLMA file does not exist: ", mlma_file)
}

# ============================================================
# 2. Trait name
# ============================================================

trait <- sub(
  "_MLMA\\.mlma$",
  "",
  basename(mlma_file)
)

cat("============================================\n")
cat("Trait:", trait, "\n")
cat("MLMA :", mlma_file, "\n")
cat("============================================\n")

# ============================================================
# 3. Output directories
# ============================================================

plot_dir <- file.path(
  base_dir,
  "plots_two_thresholds"
)

table_dir <- file.path(
  base_dir,
  "per_trait_tables"
)

dir.create(
  plot_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  table_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# 4. Thresholds
# ============================================================

n_snps <- 4549783

threshold_5e8 <- 5e-8

threshold_bonf <- 0.05 / n_snps

cat(
  "5e-8 threshold:",
  format(threshold_5e8, scientific = TRUE),
  "\n"
)

cat(
  "Bonferroni threshold:",
  format(threshold_bonf, scientific = TRUE),
  "\n"
)

cat(
  "-log10 Bonferroni threshold:",
  -log10(threshold_bonf),
  "\n"
)

# ============================================================
# 5. Read MLMA file
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
    paste(colnames(df), collapse = ", ")
  )
}

# Remove invalid P values
df <- df[
  !is.na(p) &
  is.finite(p) &
  p > 0 &
  p <= 1
]

if (nrow(df) == 0) {
  stop("No valid SNPs found for trait: ", trait)
}

cat("Valid SNPs:", nrow(df), "\n")
cat("Minimum P:", format(min(df$p), scientific = TRUE), "\n")

# ============================================================
# 6. Prepare plotting data
# ============================================================

gwas <- data.frame(
  SNP = df$SNP,
  CHR = suppressWarnings(as.integer(df$Chr)),
  BP  = df$bp,
  P   = df$p
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
# 7. Manhattan plot: P < 5e-8
# ============================================================

old_wd <- getwd()
setwd(plot_dir)

CMplot(
  gwas,
  plot.type = "m",
  LOG10 = TRUE,
  cex = 0.4,

  threshold = threshold_5e8,
  threshold.col = "red",

  chr.den.col = NULL,
  signal.col = "red",

  file = "png",
  file.output = TRUE,
  dpi = 300,

  file.name = paste0(
    trait,
    "_threshold_5e-8"
  ),

  width = 14,
  height = 6
)

setwd(old_wd)

# ============================================================
# 8. Manhattan plot: Bonferroni threshold
# ============================================================

setwd(plot_dir)

CMplot(
  gwas,
  plot.type = "m",
  LOG10 = TRUE,
  cex = 0.4,

  threshold = threshold_bonf,
  threshold.col = "red",

  chr.den.col = NULL,
  signal.col = "red",

  file = "png",
  file.output = TRUE,
  dpi = 300,

  file.name = paste0(
    trait,
    "_threshold_Bonferroni"
  ),

  width = 14,
  height = 6
)

setwd(old_wd)

# ============================================================
# 9. QQ plot
#
# No qqman package required.
# Standard expected vs observed -log10(P)
# ============================================================

qq_file <- file.path(
  plot_dir,
  paste0(
    trait,
    "_QQ.png"
  )
)

p_values <- sort(
  gwas$P
)

n <- length(
  p_values
)

expected_p <- ppoints(
  n
)

expected <- -log10(
  expected_p
)

observed <- -log10(
  p_values
)

png(
  qq_file,
  width = 1800,
  height = 1800,
  res = 220
)

plot(
  expected,
  observed,
  pch = 20,
  cex = 0.7,
  xlab = expression(
    Expected ~ -log[10](italic(P))
  ),
  ylab = expression(
    Observed ~ -log[10](italic(P))
  ),
  main = paste0(
    trait,
    " GWAS - QQ plot"
  )
)

abline(
  0,
  1,
  lty = 2
)

dev.off()

# ============================================================
# 10. Significant SNP table
#
# Since Bonferroni is stricter than 5e-8,
# selecting P < 5e-8 captures SNPs significant
# under either threshold.
# ============================================================

sig <- copy(
  df[p < threshold_5e8]
)

sig[, Trait := trait]

sig[, Significant_5e_8 :=
      p < threshold_5e8]

sig[, Significant_Bonferroni :=
      p < threshold_bonf]

sig[, Threshold_category :=
      ifelse(
        Significant_Bonferroni,
        "Bonferroni_and_5e-8",
        "5e-8_only"
      )]

sig[, Threshold_5e_8 :=
      threshold_5e8]

sig[, Threshold_Bonferroni :=
      threshold_bonf]

# Put important columns first

first_cols <- c(
  "Trait",
  "Chr",
  "SNP",
  "bp",
  "A1",
  "A2",
  "freq",
  "b",
  "se",
  "p",
  "Significant_5e_8",
  "Significant_Bonferroni",
  "Threshold_category",
  "Threshold_5e_8",
  "Threshold_Bonferroni"
)

first_cols <- first_cols[
  first_cols %in% colnames(sig)
]

remaining_cols <- setdiff(
  colnames(sig),
  first_cols
)

setcolorder(
  sig,
  c(
    first_cols,
    remaining_cols
  )
)

sig_file <- file.path(
  table_dir,
  paste0(
    trait,
    "_significant_SNPs.txt"
  )
)

fwrite(
  sig,
  sig_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

# ============================================================
# 11. Trait summary
# ============================================================

top_i <- which.min(
  df$p
)

summary_dt <- data.table(

  Trait = trait,

  Total_tested_SNPs =
    nrow(df),

  Top_SNP =
    df$SNP[top_i],

  Top_CHR =
    df$Chr[top_i],

  Top_BP =
    df$bp[top_i],

  Top_P =
    df$p[top_i],

  N_significant_5e_8 =
    sum(
      df$p < threshold_5e8
    ),

  N_significant_Bonferroni =
    sum(
      df$p < threshold_bonf
    )
)

summary_file <- file.path(
  table_dir,
  paste0(
    trait,
    "_summary.txt"
  )
)

fwrite(
  summary_dt,
  summary_file,
  sep = "\t",
  quote = FALSE
)

# ============================================================
# 12. Finish
# ============================================================

cat("\nFinished:", trait, "\n")

cat(
  "Top SNP:",
  summary_dt$Top_SNP,
  "\n"
)

cat(
  "Top P:",
  format(
    summary_dt$Top_P,
    scientific = TRUE
  ),
  "\n"
)

cat(
  "SNPs P < 5e-8:",
  summary_dt$N_significant_5e_8,
  "\n"
)

cat(
  "SNPs Bonferroni:",
  summary_dt$N_significant_Bonferroni,
  "\n"
)

cat("============================================\n")