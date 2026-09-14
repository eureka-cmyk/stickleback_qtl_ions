# =========================
# Manhattan + QQ + summary table (EDA25 covar GWAS)
# =========================

library(data.table)
library(CMplot)
library(qqman)

# -------------------------
# 1) Base directory 
# -------------------------
base_dir <- "/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_EDA10"

# dirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
dirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
dirs <- dirs[grepl("(\\/|\\\\)(As|Mg|Co|Cd|Sr|Pb|K|Fe|Ca|P|V)(_log|_no_outlier)?\\/?$", dirs)]
cat("Found", length(dirs), "GWAS trait directories under:\n", base_dir, "\n")

# -------------------------
# 2) thresholds + outputs
# -------------------------
sig_threshold  <- 5e-8
sugg_threshold <- 1e-6

out_summary <- file.path(base_dir, "All_EDA10_TopSNPs.txt")

summary_rows <- list()

# -------------------------
# 3) Loop: plot + collect hits
# -------------------------
for (d in dirs) {
  
  trait <- basename(d)
  cat("\n==============================\n")
  cat("Processing trait folder:", trait, "\n")
  cat("==============================\n")
  
  mlma_files <- list.files(d, pattern = "_MLMA\\.mlma$", full.names = TRUE)
  
  if (length(mlma_files) == 0) {
    cat("  No *_MLMA.mlma found, skip.\n")
    summary_rows[[trait]] <- data.table(
      Trait = trait, SNP = NA_character_, CHR = NA_integer_, BP = NA_integer_,
      P = NA_real_, Category = "missing_MLMA"
    )
    next
  }
  
  mlma_file <- mlma_files[1]
  cat("  Using MLMA file:", mlma_file, "\n")
  
  df <- fread(mlma_file)
  
  
  gwas <- data.frame(
    SNP = df$SNP,
    CHR = as.integer(df$Chr),
    BP  = as.integer(df$bp),
    P   = as.numeric(df$p)
  )
  gwas <- gwas[!is.na(gwas$CHR) & !is.na(gwas$BP) & !is.na(gwas$P), ]
  cat("  Loaded SNP records:", nrow(gwas), "\n")
  
  if (nrow(gwas) == 0) {
    cat("  Empty after filtering NA, skip plots.\n")
    summary_rows[[trait]] <- data.table(
      Trait = trait, SNP = NA_character_, CHR = NA_integer_, BP = NA_integer_,
      P = NA_real_, Category = "empty_after_filter"
    )
    next
  }
  
  # -------------------------
  # # 3a) Manhattan plot (CMplot) 
  # # -------------------------
  # old_wd <- getwd()
  # setwd(d)
  # 
  # CMplot(
  #   gwas,
  #   plot.type = "m",
  #   LOG10 = TRUE,
  #   cex = 0.4,
  #   threshold = c(sig_threshold, sugg_threshold),
  #   threshold.col = c("red", "blue"),
  #   chr.den.col = NULL,
  #   signal.col = "red",
  #   file = "png",
  #   file.output = TRUE,
  #   dpi = 300,
  #   file.name = paste0(trait, "_EDA25"),
  #   width = 14,
  #   height = 6
  # )
  # 
  # # -------------------------
  # # 3b) QQ plot (qqman)
  # # -------------------------
  # out_qq <- file.path(d, paste0(trait, "_EDA25_QQ.png"))
  # png(out_qq, width = 1600, height = 1600, res = 200)
  # qq(gwas$P, main = paste("QQ Plot —", trait, "(EDA25 covar)"))
  # dev.off()
  # 
  # setwd(old_wd)
  # 
  # cat("  Plots saved in:", d, "\n")
  # 
  # -------------------------
  # 3c) Collect top SNPs for summary
  # -------------------------
  sig_snps  <- gwas[gwas$P < sig_threshold, ]
  sugg_snps <- gwas[gwas$P < sugg_threshold & gwas$P >= sig_threshold, ]
  
  if (nrow(sig_snps) > 0) {
    summary_rows[[paste0(trait, "_sig")]] <- data.table(
      Trait = trait,
      SNP = sig_snps$SNP,
      CHR = sig_snps$CHR,
      BP  = sig_snps$BP,
      P   = sig_snps$P,
      Category = "significant"
    )
  }
  
  if (nrow(sugg_snps) > 0) {
    summary_rows[[paste0(trait, "_sugg")]] <- data.table(
      Trait = trait,
      SNP = sugg_snps$SNP,
      CHR = sugg_snps$CHR,
      BP  = sugg_snps$BP,
      P   = sugg_snps$P,
      Category = "suggestive"
    )
  }
  
  if (nrow(sig_snps) == 0 && nrow(sugg_snps) == 0) {
    summary_rows[[paste0(trait, "_nohit")]] <- data.table(
      Trait = trait, SNP = NA_character_, CHR = NA_integer_, BP = NA_integer_,
      P = NA_real_, Category = "no_hit"
    )
  }
}




# -------------------------
# 4) Write summary table
# -------------------------
summary_dt <- rbindlist(summary_rows, use.names = TRUE, fill = TRUE)
setorder(summary_dt, Trait, Category, P)

fwrite(summary_dt, out_summary, sep = "\t", quote = FALSE, na = "NA")
cat("\nAll top SNPs saved to:\n", out_summary, "\n")


# =========================
# Combine Manhattan PNGs -> PDF
# + Split summary table into (11 elements) and (6 main)
# 解决：Trait 名含 _log/_no_outlier 导致筛选为空
# =========================

library(data.table)
library(grid)
library(png)

base_dir <- "/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_EDA10"


# -------------------------
clean_element <- function(x) {
  x <- sub("_no_outlier$", "", x)
  x <- sub("_log$", "", x)
  x
}


summary_dt[, Element := clean_element(Trait)]

eleven_elements <- c("As","Mg","Co","Cd","Sr","Pb","K","Fe","Ca","P","V")
six_elements    <- c("Ca","Co","Fe","Sr","P","V")

out_eleven_table <- file.path(base_dir, "nine_log_no_outlier_eda.txt")  # 按你命名，但内容含P/V共11个
out_six_table    <- file.path(base_dir, "six_main_eda.txt")

eleven_table <- summary_dt[Element %in% eleven_elements]
six_table    <- summary_dt[Element %in% six_elements]

setorder(eleven_table, Element, Category, P)
setorder(six_table,    Element, Category, P)

fwrite(eleven_table, out_eleven_table, sep = "\t", quote = FALSE, na = "NA")
fwrite(six_table,    out_six_table,    sep = "\t", quote = FALSE, na = "NA")

cat("Saved 11-element table (named nine_log_no_outlier_10kbeda):\n", out_eleven_table, "\n")
cat("Saved 6-element table:\n", out_six_table, "\n")


trait_dirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
trait_names <- basename(trait_dirs)

dir_dt <- data.table(
  TraitFolder = trait_names,
  Dir = trait_dirs
)
dir_dt[, Element := clean_element(TraitFolder)]

find_manhattan_png_in_dir <- function(d) {
  pngs <- list.files(d, pattern = "\\.png$", full.names = TRUE)
  pngs <- pngs[!grepl("_QQ\\.png$", pngs)]
  eda <- pngs[grepl("EDA+-10kb", basename(pngs), ignore.case = TRUE)]
  if (length(eda) > 0) return(eda[1])
  if (length(pngs) > 0) return(pngs[1])
  NA_character_
}


pick_dir_for_element <- function(el) {
  cand <- dir_dt[Element == el]
  if (nrow(cand) == 0) return(NULL)
  cand_no <- cand[grepl("_no_outlier$", TraitFolder)]
  if (nrow(cand_no) > 0) return(cand_no[1])
  cand_log <- cand[grepl("_log$", TraitFolder)]
  if (nrow(cand_log) > 0) return(cand_log[1])
  cand[1]
}

get_ordered_pngs <- function(elements) {
  out <- data.table(Element = elements, TraitFolder = NA_character_, Dir = NA_character_, Png = NA_character_)
  for (i in seq_along(elements)) {
    el <- elements[i]
    pick <- pick_dir_for_element(el)
    if (!is.null(pick)) {
      out[i, TraitFolder := pick$TraitFolder]
      out[i, Dir := pick$Dir]
      out[i, Png := find_manhattan_png_in_dir(pick$Dir)]
    }
  }
  out
}

draw_panel <- function(img_path, label) {
  pushViewport(viewport(layout = grid.layout(2, 1, heights = unit(c(1.2, 8.8), "null"))))
  
  pushViewport(viewport(layout.pos.row = 1))
  grid.text(label, gp = gpar(cex = 1.8, fontface = "bold")) 
  popViewport()
  
  pushViewport(viewport(layout.pos.row = 2))
  if (!is.na(img_path) && file.exists(img_path)) {
    img <- readPNG(img_path)
    grid.raster(img, interpolate = TRUE) 
  } else {
    grid.rect(gp = gpar(col = "grey70", lty = "dashed"))
    grid.text(paste0("(Missing: ", label, ")"), gp = gpar(cex = 1, col = "grey50"))
  }
  popViewport()
  
  popViewport()
}

make_pdf_multi_page <- function(png_dt, outpdf, ncol = 2, nrow = 4, width = 20, height = 16) {
  n <- nrow(png_dt)
  plots_per_page <- ncol * nrow
  
  pdf(outpdf, width = width, height = height)
  
  for (i in seq_len(n)) {
    if ((i - 1) %% plots_per_page == 0) {
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(nrow, ncol)))
    }
    
    idx_in_page <- (i - 1) %% plots_per_page
    curr_r <- (idx_in_page %/% ncol) + 1
    curr_c <- (idx_in_page %% ncol) + 1
    
    pushViewport(viewport(layout.pos.row = curr_r, layout.pos.col = curr_c))
    pushViewport(viewport(width = 0.95, height = 0.95)) 
    draw_panel(png_dt$Png[i], png_dt$Element[i])
    popViewport(2)
    
    if (i %% plots_per_page == 0 || i == n) {
      popViewport()
    }
  }
  
  dev.off()
  cat("Saved formatted PDF:\n", outpdf, "\n")
}

eleven_pngs <- get_ordered_pngs(eleven_elements)
outpdf_eleven <- file.path(base_dir, "nine_log_no_outlier_10kbeda_manhattan.pdf")
make_pdf_multi_page(eleven_pngs, outpdf_eleven, ncol = 2, nrow = 3, width = 18, height = 16)

six_pngs <- get_ordered_pngs(six_elements)
outpdf_six <- file.path(base_dir, "six_main_eda_manhattan.pdf")
make_pdf_multi_page(six_pngs, outpdf_six, ncol = 2, nrow = 3, width = 18, height = 16)

cat("\n[Check missing Manhattan PNGs]\n")
print(eleven_pngs[is.na(Png) | !file.exists(Png), .(Element, TraitFolder, Dir, Png)])
print(six_pngs[is.na(Png) | !file.exists(Png), .(Element, TraitFolder, Dir, Png)])

