#!/usr/bin/env Rscript

# ============================================================
# Merge no-PC GWAS result tables
#
# Input:
#   per_trait_tables/*_summary.txt
#   per_trait_tables/*_significant_SNPs.txt
#
# Output:
#   1. Significant_SNPs_two_thresholds.txt
#   2. Significant_SNP_counts_two_thresholds.txt
#   3. GWAS_significance_thresholds.txt
#
# No plotting or PDF generation in this script.
# ============================================================

.libPaths(c(
  "/gpfs01/home/shxyt1/R/x86_64-pc-linux-gnu-library/4.4",
  .libPaths()
))

suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# 1. Paths
# ============================================================

base_dir <- "/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_h_w_log_no_outlier"

table_dir <- file.path(
  base_dir,
  "per_trait_tables"
)

if (!dir.exists(table_dir)) {
  stop(
    "Directory does not exist: ",
    table_dir
  )
}

# ============================================================
# 2. Thresholds
# ============================================================

n_snps <- 4549783

threshold_5e8 <- 5e-8

threshold_bonf <- 0.05 / n_snps

cat("============================================\n")
cat("Number of SNPs:", n_snps, "\n")
cat(
  "P < 5e-8 threshold:",
  format(threshold_5e8, scientific = TRUE),
  "\n"
)
cat(
  "Bonferroni threshold:",
  format(threshold_bonf, scientific = TRUE),
  "\n"
)
cat("============================================\n")


# ============================================================
# 3. Merge trait summary tables
# ============================================================

summary_files <- list.files(
  table_dir,
  pattern = "_summary\\.txt$",
  full.names = TRUE
)

cat(
  "\nSummary files found:",
  length(summary_files),
  "\n"
)

if (length(summary_files) > 0) {

  summary_list <- lapply(
    summary_files,
    fread
  )

  summary_table <- rbindlist(
    summary_list,
    fill = TRUE,
    use.names = TRUE
  )

  if ("Trait" %in% names(summary_table)) {
    setorder(
      summary_table,
      Trait
    )
  }

  summary_out <- file.path(
    base_dir,
    "Significant_SNP_counts_two_thresholds.txt"
  )

  fwrite(
    summary_table,
    summary_out,
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )

  cat(
    "Combined summary table written:\n",
    summary_out,
    "\n"
  )

} else {

  warning(
    "No *_summary.txt files found."
  )
}


# ============================================================
# 4. Merge significant SNP tables
# ============================================================

sig_files <- list.files(
  table_dir,
  pattern = "_significant_SNPs\\.txt$",
  full.names = TRUE
)

cat(
  "\nSignificant SNP files found:",
  length(sig_files),
  "\n"
)

sig_list <- list()

for (f in sig_files) {

  tmp <- fread(
    f,
    header = TRUE
  )

  # Ignore traits with zero significant SNPs
  if (nrow(tmp) > 0) {

    sig_list[[
      length(sig_list) + 1
    ]] <- tmp
  }
}


sig_out <- file.path(
  base_dir,
  "Significant_SNPs_two_thresholds.txt"
)


if (length(sig_list) > 0) {

  sig_table <- rbindlist(
    sig_list,
    fill = TRUE,
    use.names = TRUE
  )

  # Sort by element and P value
  if (
    all(
      c("Trait", "p") %in%
        names(sig_table)
    )
  ) {

    setorder(
      sig_table,
      Trait,
      p
    )
  }

  # Put important columns first
  preferred_cols <- c(
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

  preferred_cols <- preferred_cols[
    preferred_cols %in%
      names(sig_table)
  ]

  other_cols <- setdiff(
    names(sig_table),
    preferred_cols
  )

  setcolorder(
    sig_table,
    c(
      preferred_cols,
      other_cols
    )
  )

  fwrite(
    sig_table,
    sig_out,
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )

  cat(
    "\nCombined significant SNP table written:\n",
    sig_out,
    "\n"
  )

  cat(
    "Total significant SNP rows:",
    nrow(sig_table),
    "\n"
  )

  # ----------------------------------------------------------
  # Print counts for each threshold
  # ----------------------------------------------------------

  if (
    "Significant_5e_8" %in%
      names(sig_table)
  ) {

    cat(
      "SNP rows significant at P < 5e-8:",
      sum(
        sig_table$Significant_5e_8,
        na.rm = TRUE
      ),
      "\n"
    )
  }

  if (
    "Significant_Bonferroni" %in%
      names(sig_table)
  ) {

    cat(
      "SNP rows significant at Bonferroni threshold:",
      sum(
        sig_table$Significant_Bonferroni,
        na.rm = TRUE
      ),
      "\n"
    )
  }

} else {

  cat(
    "\nNo significant SNPs were found in the per-trait tables.\n"
  )
}


# ============================================================
# 5. Write threshold information
# ============================================================

threshold_table <- data.table(

  Method = c(
    "Standard_5e-8",
    "Bonferroni"
  ),

  Number_of_SNPs = c(
    n_snps,
    n_snps
  ),

  P_threshold = c(
    threshold_5e8,
    threshold_bonf
  ),

  minus_log10_P = c(
    -log10(threshold_5e8),
    -log10(threshold_bonf)
  )
)

threshold_out <- file.path(
  base_dir,
  "GWAS_significance_thresholds.txt"
)

fwrite(
  threshold_table,
  threshold_out,
  sep = "\t",
  quote = FALSE
)

cat(
  "\nThreshold information written:\n",
  threshold_out,
  "\n"
)


# ============================================================
# 6. Finish
# ============================================================

cat("\n")
cat("============================================\n")
cat("TABLE MERGE FINISHED\n")
cat("============================================\n")

cat(
  "Bonferroni threshold =",
  format(
    threshold_bonf,
    digits = 7,
    scientific = TRUE
  ),
  "\n"
)