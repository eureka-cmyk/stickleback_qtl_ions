#!/usr/bin/env Rscript

.libPaths(c(
  "/gpfs01/home/shxyt1/R/x86_64-pc-linux-gnu-library/4.4",
  .libPaths()
))

suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# 1. Thresholds
# ============================================================

n_snps <- 4549783

threshold_5e8  <- 5e-8
threshold_bonf <- 0.05 / n_snps

cat("5e-8 threshold:", threshold_5e8, "\n")
cat("Bonferroni threshold:", threshold_bonf, "\n")


# ============================================================
# 2. Significant SNP summary tables
# ============================================================

sig_5pc_file <- paste0(
  "/gpfs01/home/shxyt1/data/sticklebacks/gwas/",
  "gwas_5pc_log_nine_main/5PC_significant_summary/",
  "Significant_SNPs_two_thresholds_5PC.txt"
)

sig_nopc_file <- paste0(
  "/gpfs01/home/shxyt1/data/sticklebacks/gwas/",
  "gwas_h_w_log_no_outlier/",
  "Significant_SNPs_two_thresholds.txt"
)


# ============================================================
# 3. Full MLMA files
# ============================================================

mlma_5pc <- c(

  Ca = paste0(
    "/gpfs01/home/shxyt1/data/sticklebacks/gwas/",
    "gwas_5pc_log_nine_main/Ca_no_outlier/",
    "Ca_no_outlier_MLMA.mlma"
  ),

  Sr = paste0(
    "/gpfs01/home/shxyt1/data/sticklebacks/gwas/",
    "gwas_5pc_log_nine_main/Sr_no_outlier/",
    "Sr_no_outlier_MLMA.mlma"
  ),

  Co = paste0(
    "/gpfs01/home/shxyt1/data/sticklebacks/gwas/",
    "gwas_5pc_log_nine_main/Co_no_outlier/",
    "Co_no_outlier_MLMA.mlma"
  ),

  Fe = paste0(
    "/gpfs01/home/shxyt1/data/sticklebacks/gwas/",
    "gwas_5pc_log_nine_main/Fe_no_outlier/",
    "Fe_no_outlier_MLMA.mlma"
  ),

  P = paste0(
    "/gpfs01/home/shxyt1/data/sticklebacks/gwas/",
    "gwas_5pc_log_no_outlier/P_log/",
    "P_log_MLMA.mlma"
  ),

  V = paste0(
    "/gpfs01/home/shxyt1/data/sticklebacks/gwas/",
    "gwas_5pc_log_no_outlier/V_log/",
    "V_log_MLMA.mlma"
  )
)


mlma_nopc <- c(

  Ca = paste0(
    "/gpfs01/home/shxyt1/data/sticklebacks/gwas/",
    "gwas_h_w_log_no_outlier/Ca/",
    "Ca_MLMA.mlma"
  ),

  Sr = paste0(
    "/gpfs01/home/shxyt1/data/sticklebacks/gwas/",
    "gwas_h_w_log_no_outlier/Sr/",
    "Sr_MLMA.mlma"
  ),

  Co = paste0(
    "/gpfs01/home/shxyt1/data/sticklebacks/gwas/",
    "gwas_h_w_log_no_outlier/Co/",
    "Co_MLMA.mlma"
  ),

  Fe = paste0(
    "/gpfs01/home/shxyt1/data/sticklebacks/gwas/",
    "gwas_h_w_log_no_outlier/Fe/",
    "Fe_MLMA.mlma"
  ),

  P = paste0(
    "/gpfs01/home/shxyt1/data/sticklebacks/gwas/",
    "gwas_h_w_log_no_outlier/P/",
    "P_MLMA.mlma"
  ),

  V = paste0(
    "/gpfs01/home/shxyt1/data/sticklebacks/gwas/",
    "gwas_h_w_log_no_outlier/V/",
    "V_MLMA.mlma"
  )
)


traits <- c(
  "Ca",
  "Sr",
  "P",
  "Co",
  "Fe",
  "V"
)


# ============================================================
# 4. Output directory
# ============================================================

out_dir <- paste0(
  "/gpfs01/home/shxyt1/data/sticklebacks/gwas/",
  "gwas_h_w_log_no_outlier/5PC_vs_noPC"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 5. Read significant SNP tables
# ============================================================

sig5 <- fread(sig_5pc_file)
sign <- fread(sig_nopc_file)

sig5 <- sig5[Trait %in% traits]
sign <- sign[Trait %in% traits]

cat("\n5-PC significant rows:", nrow(sig5), "\n")
cat("no-PC significant rows:", nrow(sign), "\n")


# ============================================================
# 6. Union of SNPs significant in either model
# ============================================================

candidate_snps <- unique(
  rbind(
    sig5[, .(Trait, SNP)],
    sign[, .(Trait, SNP)]
  )
)

setorder(
  candidate_snps,
  Trait,
  SNP
)

cat(
  "Union of significant SNPs:",
  nrow(candidate_snps),
  "\n"
)


# ============================================================
# 7. Read selected SNPs from one MLMA file
# ============================================================

read_selected_snps <- function(
  mlma_file,
  snps,
  model
) {

  if (!file.exists(mlma_file)) {
    stop("MLMA file not found: ", mlma_file)
  }

  d <- fread(
    mlma_file,
    select = c(
      "Chr",
      "SNP",
      "bp",
      "A1",
      "A2",
      "Freq",
      "b",
      "se",
      "p"
    ),
    showProgress = FALSE
  )

  d <- d[
    SNP %in% snps
  ]

  if (model == "5PC") {

    setnames(
      d,
      old = c(
        "Chr",
        "bp",
        "A1",
        "A2",
        "Freq",
        "b",
        "se",
        "p"
      ),
      new = c(
        "Chr_5PC",
        "BP_5PC",
        "A1_5PC",
        "A2_5PC",
        "Freq_5PC",
        "Beta_5PC",
        "SE_5PC",
        "P_5PC"
      )
    )

  } else {

    setnames(
      d,
      old = c(
        "Chr",
        "bp",
        "A1",
        "A2",
        "Freq",
        "b",
        "se",
        "p"
      ),
      new = c(
        "Chr_noPC",
        "BP_noPC",
        "A1_noPC",
        "A2_noPC",
        "Freq_noPC",
        "Beta_noPC",
        "SE_noPC",
        "P_noPC"
      )
    )
  }

  return(d)
}


# ============================================================
# 8. Compare each trait
# ============================================================

comparison_list <- list()

for (trait in traits) {

  cat("\n============================================\n")
  cat("Trait:", trait, "\n")
  cat("============================================\n")

  snps <- candidate_snps[
    Trait == trait,
    SNP
  ]

  cat(
    "Candidate SNPs:",
    length(snps),
    "\n"
  )

  if (length(snps) == 0) {
    next
  }


  # ----------------------------------------------------------
  # Read same SNPs from 5-PC model
  # ----------------------------------------------------------

  cat("Reading 5-PC MLMA...\n")

  d5 <- read_selected_snps(
    mlma_5pc[[trait]],
    snps,
    "5PC"
  )

  cat(
    "Found in 5-PC:",
    nrow(d5),
    "\n"
  )


  # ----------------------------------------------------------
  # Read same SNPs from no-PC model
  # ----------------------------------------------------------

  cat("Reading no-PC MLMA...\n")

  dn <- read_selected_snps(
    mlma_nopc[[trait]],
    snps,
    "noPC"
  )

  cat(
    "Found in no-PC:",
    nrow(dn),
    "\n"
  )


  # ----------------------------------------------------------
  # Merge
  # ----------------------------------------------------------

  comp <- merge(
    d5,
    dn,
    by = "SNP",
    all = TRUE
  )

  comp[, Trait := trait]


  # ----------------------------------------------------------
  # Final chromosome / position
  # ----------------------------------------------------------

  comp[
    ,
    Chr := ifelse(
      !is.na(Chr_5PC),
      Chr_5PC,
      Chr_noPC
    )
  ]

  comp[
    ,
    BP := ifelse(
      !is.na(BP_5PC),
      BP_5PC,
      BP_noPC
    )
  ]


  # ==========================================================
  # 9. Significance flags: 5e-8
  # ==========================================================

  comp[
    ,
    Sig_5PC_5e8 :=
      !is.na(P_5PC) &
      P_5PC < threshold_5e8
  ]

  comp[
    ,
    Sig_noPC_5e8 :=
      !is.na(P_noPC) &
      P_noPC < threshold_5e8
  ]


  # ==========================================================
  # 10. Significance flags: Bonferroni
  # ==========================================================

  comp[
    ,
    Sig_5PC_Bonf :=
      !is.na(P_5PC) &
      P_5PC < threshold_bonf
  ]

  comp[
    ,
    Sig_noPC_Bonf :=
      !is.na(P_noPC) &
      P_noPC < threshold_bonf
  ]


  # ==========================================================
  # 11. Classification at 5e-8
  # ==========================================================

  comp[
    ,
    Comparison_5e8 :=
      ifelse(
        Sig_5PC_5e8 & Sig_noPC_5e8,
        "Both",
        ifelse(
          Sig_5PC_5e8 & !Sig_noPC_5e8,
          "5PC_only",
          ifelse(
            !Sig_5PC_5e8 & Sig_noPC_5e8,
            "noPC_only",
            "Neither"
          )
        )
      )
  ]


  # ==========================================================
  # 12. Classification at Bonferroni threshold
  # ==========================================================

  comp[
    ,
    Comparison_Bonf :=
      ifelse(
        Sig_5PC_Bonf & Sig_noPC_Bonf,
        "Both",
        ifelse(
          Sig_5PC_Bonf & !Sig_noPC_Bonf,
          "5PC_only",
          ifelse(
            !Sig_5PC_Bonf & Sig_noPC_Bonf,
            "noPC_only",
            "Neither"
          )
        )
      )
  ]


  # ==========================================================
  # 13. Effect-size comparison
  # ==========================================================

  comp[
    ,
    Beta_change :=
      Beta_noPC - Beta_5PC
  ]

  comp[
    ,
    Beta_change_percent :=
      100 *
      (Beta_noPC - Beta_5PC) /
      abs(Beta_5PC)
  ]

  comp[
    ,
    SE_change :=
      SE_noPC - SE_5PC
  ]

  comp[
    ,
    SE_change_percent :=
      100 *
      (SE_noPC / SE_5PC - 1)
  ]

  comp[
    ,
    P_ratio_noPC_vs_5PC :=
      P_noPC / P_5PC
  ]


  # ==========================================================
  # 14. Allele checks
  # ==========================================================

  comp[
    ,
    Same_A1 :=
      A1_5PC == A1_noPC
  ]

  comp[
    ,
    Same_A2 :=
      A2_5PC == A2_noPC
  ]


  # ==========================================================
  # 15. Save trait result
  # ==========================================================

  comparison_list[[trait]] <- comp


  rm(
    d5,
    dn,
    comp
  )

  gc()
}


# ============================================================
# 16. Combine all traits
# ============================================================

final <- rbindlist(
  comparison_list,
  fill = TRUE,
  use.names = TRUE
)


# ============================================================
# 17. Column order
# ============================================================

final_cols <- c(

  "Trait",
  "Chr",
  "BP",
  "SNP",

  "A1_5PC",
  "A2_5PC",
  "Freq_5PC",

  "A1_noPC",
  "A2_noPC",
  "Freq_noPC",

  "Beta_5PC",
  "SE_5PC",
  "P_5PC",

  "Beta_noPC",
  "SE_noPC",
  "P_noPC",

  "Beta_change",
  "Beta_change_percent",

  "SE_change",
  "SE_change_percent",

  "P_ratio_noPC_vs_5PC",

  "Sig_5PC_5e8",
  "Sig_noPC_5e8",
  "Comparison_5e8",

  "Sig_5PC_Bonf",
  "Sig_noPC_Bonf",
  "Comparison_Bonf",

  "Same_A1",
  "Same_A2"
)


final <- final[
  ,
  ..final_cols
]


setorder(
  final,
  Trait,
  Chr,
  BP
)


# ============================================================
# 18. Write detailed comparison table
# ============================================================

out_file <- file.path(
  out_dir,
  "GWAS_5PC_vs_noPC_significant_SNP_comparison.txt"
)

fwrite(
  final,
  out_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ============================================================
# 19. Trait-level summary
# ============================================================

summary_table <- final[
  ,
  .(

    N_union_SNPs = .N,

    N_both_5e8 =
      sum(
        Comparison_5e8 == "Both",
        na.rm = TRUE
      ),

    N_5PC_only_5e8 =
      sum(
        Comparison_5e8 == "5PC_only",
        na.rm = TRUE
      ),

    N_noPC_only_5e8 =
      sum(
        Comparison_5e8 == "noPC_only",
        na.rm = TRUE
      ),

    N_both_Bonf =
      sum(
        Comparison_Bonf == "Both",
        na.rm = TRUE
      ),

    N_5PC_only_Bonf =
      sum(
        Comparison_Bonf == "5PC_only",
        na.rm = TRUE
      ),

    N_noPC_only_Bonf =
      sum(
        Comparison_Bonf == "noPC_only",
        na.rm = TRUE
      )

  ),
  by = Trait
]


setorder(
  summary_table,
  Trait
)


summary_file <- file.path(
  out_dir,
  "GWAS_5PC_vs_noPC_trait_summary.txt"
)

fwrite(
  summary_table,
  summary_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)


# ============================================================
# 20. Check Ca chr14:321923
# ============================================================

cat("\n============================================\n")
cat("Check: Ca 14:321923\n")
cat("============================================\n")

print(
  final[
    Trait == "Ca" &
    SNP == "14:321923"
  ]
)


# ============================================================
# 21. Print summary
# ============================================================

cat("\n============================================\n")
cat("Trait summary\n")
cat("============================================\n")

print(
  summary_table
)


cat("\n============================================\n")
cat("FINISHED\n")
cat("============================================\n")

cat("\nDetailed table:\n")
cat(out_file, "\n")

cat("\nSummary table:\n")
cat(summary_file, "\n")