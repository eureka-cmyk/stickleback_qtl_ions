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
# Load GWAS results
# ===========================
df <- fread(mlma_file)

gwas <- data.frame(
  SNP = df$SNP,
  CHR = as.integer(df$Chr),
  BP  = df$bp,
  P   = df$p
)

gwas <- gwas[complete.cases(gwas), ]
cat("Loaded SNPs:", nrow(gwas), "\n")

# ===========================
# Manhattan plot
# ===========================
CMplot(
  gwas,
  plot.type = "m",
  LOG10 = TRUE,
  cex = 0.4,
  threshold = c(5e-8, 1e-6),
  threshold.col = c("red", "blue"),
  chr.den.col = NULL,
  signal.col = "red",
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
png(
  filename = paste0(element, "_QQ.png"),
  width = 1600,
  height = 1600,
  res = 200
)
qq(gwas$P, main = paste("QQ plot —", element))
dev.off()

cat("Finished:", element, "\n")
