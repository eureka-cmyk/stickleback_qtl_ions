library(data.table)

locus_file <- "/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_nine_main/six_main/regional plot/locus_summary_detailed_stats.csv"

ld_dir <- "/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_nine_main/six_main"

plot_flank_bp <- 150000L

locus_tbl <- fread(locus_file)

locus_tbl[, CHR := as.integer(CHR)]
locus_tbl[, locus_id := as.integer(locus_id)]
locus_tbl[, lead_BP := as.integer(lead_BP)]
locus_tbl[, Element := as.character(Element)]
locus_tbl[, lead_SNP := as.character(lead_SNP)]

elements_order <- c("Ca", "Co", "Fe", "Sr", "P", "V")
locus_tbl <- locus_tbl[Element %in% elements_order]
locus_tbl <- locus_tbl[order(factor(Element, levels = elements_order), CHR, locus_id)]

check_one_locus <- function(trait, chr_target, locus_id, lead_snp, lead_bp) {
  
  start_target <- max(1L, lead_bp - plot_flank_bp)
  end_target   <- lead_bp + plot_flank_bp
  
  ld_file1 <- file.path(ld_dir, sprintf("%s_chr%d_locus%d_lead%s.ld",
                                        trait, chr_target, locus_id, lead_snp))
  
  ld_file2 <- file.path(ld_dir, sprintf("%s_chr%d_locus%d.ld",
                                        trait, chr_target, locus_id))
  
  ld_file3 <- Sys.glob(file.path(ld_dir, sprintf("%s_chr%d_locus%d*.ld",
                                                 trait, chr_target, locus_id)))
  
  if (file.exists(ld_file1)) {
    ld_file <- ld_file1
  } else if (file.exists(ld_file2)) {
    ld_file <- ld_file2
  } else if (length(ld_file3) > 0) {
    ld_file <- ld_file3[1]
  } else {
    ld_file <- NA_character_
  }
  
  cat("\n========================================\n")
  cat("Trait:", trait, "\n")
  cat("CHR:", chr_target, "\n")
  cat("locus_id:", locus_id, "\n")
  cat("lead_snp:", lead_snp, "\n")
  cat("lead_bp:", lead_bp, "\n")
  cat("plot window:", start_target, "-", end_target, "\n")
  
  cat("ld_file1:", ld_file1, "\n")
  cat("ld_file2:", ld_file2, "\n")
  cat("ld_file3 matched:", paste(ld_file3, collapse = " | "), "\n")
  cat("selected ld_file:", ld_file, "\n")
  
  if (is.na(ld_file) || !file.exists(ld_file)) {
    cat("RESULT: LD file not found.\n")
    return(NULL)
  }
  
  ld <- fread(ld_file)
  
  cat("LD columns:", paste(names(ld), collapse = ", "), "\n")
  cat("LD rows:", nrow(ld), "\n")
  
  if (!all(c("SNP_A", "SNP_B", "R2") %in% names(ld))) {
    cat("RESULT: LD file does not contain SNP_A, SNP_B, R2.\n")
    return(NULL)
  }
  
  n_lead_A <- sum(ld$SNP_A == lead_snp)
  n_lead_B <- sum(ld$SNP_B == lead_snp)
  n_lead_any <- sum(ld$SNP_A == lead_snp | ld$SNP_B == lead_snp)
  
  cat("lead in SNP_A:", n_lead_A, "\n")
  cat("lead in SNP_B:", n_lead_B, "\n")
  cat("lead in either SNP_A/SNP_B:", n_lead_any, "\n")
  
  all_snps <- unique(c(as.character(ld$SNP_A), as.character(ld$SNP_B)))
  
  snp_bp <- suppressWarnings(as.integer(sub(".*:", "", all_snps)))
  snp_bp <- snp_bp[!is.na(snp_bp)]
  
  if (length(snp_bp) > 0) {
    cat("LD SNP BP range:", min(snp_bp), "-", max(snp_bp), "\n")
    cat("LD covers plot start?", min(snp_bp) <= start_target, "\n")
    cat("LD covers plot end?", max(snp_bp) >= end_target, "\n")
  } else {
    cat("Could not parse BP from SNP names.\n")
  }
  
  if (n_lead_any == 0) {
    cat("RESULT: Lead SNP is not in this LD file. This locus will have no LD colour.\n")
  } else if (n_lead_A == 0 && n_lead_B > 0) {
    cat("RESULT: Lead SNP is only in SNP_B. Your plotting code must read SNP_B == lead_snp as well.\n")
  } else if (length(snp_bp) > 0 && (min(snp_bp) > start_target || max(snp_bp) < end_target)) {
    cat("RESULT: LD file exists, but LD range is smaller than the 150 kb plot window. Edges will be grey.\n")
  } else {
    cat("RESULT: LD file should work for this plot window.\n")
  }
}

for (i in seq_len(nrow(locus_tbl))) {
  check_one_locus(
    trait = locus_tbl$Element[i],
    chr_target = locus_tbl$CHR[i],
    locus_id = locus_tbl$locus_id[i],
    lead_snp = locus_tbl$lead_SNP[i],
    lead_bp = locus_tbl$lead_BP[i]
  )
}