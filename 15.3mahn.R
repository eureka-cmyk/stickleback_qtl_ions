library(data.table)
library(CMplot)
library(qqman)

options(bitmapType = "cairo")

# ===========================
# Parse command line arguments
# ===========================
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag) {
  idx <- which(args == flag)
  if (length(idx) == 0 || idx == length(args)) return(NA)
  args[idx + 1]
}

base_dir <- get_arg("--base")
element  <- get_arg("--element")
out_dir  <- get_arg("--outdir")

if (any(is.na(c(base_dir, element, out_dir)))) {
  stop(
    "Usage:\n",
    "Rscript manhattan.R ",
    "--base <GWAS_BASE_DIR> ",
    "--element <ELEMENT_NAME> ",
    "--outdir <OUTPUT_DIR>\n"
  )
}

# ===========================
# Paths
# ===========================
work_dir  <- file.path(base_dir, element)
mlma_file <- file.path(work_dir, paste0(element, "_MLMA.mlma"))

cat("Processing element:", element, "\n")
cat("MLMA file:", mlma_file, "\n")
cat("Output dir:", out_dir, "\n")

if (!file.exists(mlma_file)) {
  stop("MLMA file not found: ", mlma_file)
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
setwd(out_dir)

# ===========================
# Load GWAS results (GCTA MLMA)
# ===========================
df <- fread(mlma_file)

# Be robust to column naming differences across GCTA versions
col_lower <- tolower(names(df))

get_col <- function(possible) {
  # Return the first matching column name (case-insensitive)
  hit <- which(col_lower %in% tolower(possible))
  if (length(hit) == 0) return(NULL)
  names(df)[hit[1]]
}

col_snp <- get_col(c("snp"))
col_chr <- get_col(c("chr", "chrom", "chromosome"))
col_bp  <- get_col(c("bp", "pos", "position"))
col_p   <- get_col(c("p", "pvalue", "pval"))

if (is.null(col_snp) || is.null(col_chr) || is.null(col_bp) || is.null(col_p)) {
  stop("Cannot find required columns in MLMA file. Found columns: ",
       paste(names(df), collapse = ", "))
}

gwas <- data.frame(
  SNP = df[[col_snp]],
  CHR = as.integer(df[[col_chr]]),
  BP  = as.integer(df[[col_bp]]),
  P   = as.numeric(df[[col_p]])
)

gwas <- gwas[complete.cases(gwas), ]
cat("Loaded SNPs:", nrow(gwas), "\n")

if (nrow(gwas) == 0) stop("No valid GWAS rows after filtering complete cases.")

# ===========================
# Manhattan plot
# ===========================
# CMplot(
#   gwas,
#   plot.type = "m",
#   LOG10 = TRUE,
#   cex = 0.4,
#   threshold = c(5e-8, 1e-6),
#   threshold.col = c("red", "blue"),
#   chr.den.col = NULL,
#   signal.col = "red",
#   file = "png",
#   file.output = TRUE,
#   dpi = 300,
#   file.name = paste0(element),
#   width = 14,
#   height = 6
# )

# ===========================
# Manhattan plot
# ===========================
CMplot(
  gwas,
  plot.type = "m",
  LOG10 = TRUE,
  cex = 0.4,

  col = c("grey20", "grey70"),

  threshold = c(5e-8),
  threshold.col = c("red"),
  threshold.lty = c(2, 2),
  threshold.lwd = c(1.2, 1.2),
  
  chr.den.col = NULL,
  signal.col = "red",
  signal.cex = 0.6,
  
  file = "png",
  file.output = TRUE,
  dpi = 300,
  file.name = paste0(element),
  width = 14,
  height = 6
)

# ===========================
# QQ plot
# ===========================
# png(
#   filename = paste0(element, "_QQ.png"),
#   width = 1600,
#   height = 1600,
#   res = 200
# )
# qq(gwas$P, main = paste("QQ plot —", element))
dev.off()

cat("Finished:", element, "\n")
