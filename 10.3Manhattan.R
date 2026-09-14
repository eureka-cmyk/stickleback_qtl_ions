# ===========================
#   Load libraries
# ===========================
library(data.table)
library(CMplot)
library(qqman)

# ===========================
#   File paths
# ===========================
mlma_file <- "/gpfs01/home/shxyt1/data/sticklebacks/gwas/height_weight/height/height_MLMA.mlma"

outdir <- "/gpfs01/home/shxyt1/data/sticklebacks/gwas/height_weight/height"
setwd(outdir)
out_qq <- file.path(outdir, "height_QQ.png")

# ===========================
#   Load MLMA results
# ===========================
df <- fread(mlma_file, header = TRUE)

head(df)
# Make sure columns exist
if(!("p" %in% names(df))) stop("Cannot find column 'p' in MLMA file!")
if(!("Chr" %in% names(df))) stop("Cannot find column 'Chr'!")

# Extract columns in PLINK-like format
gwas <- data.frame(
  SNP = df$SNP,
  CHR = as.integer(df$Chr),
  BP  = df$bp,
  P   = df$p
)

# Drop SNPs with invalid data
gwas <- gwas[!is.na(gwas$CHR) & !is.na(gwas$BP) & !is.na(gwas$P), ]

cat("Loaded GWAS records:", nrow(gwas), "\n")

# ===========================
#   Manhattan plot (CMplot)
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
#  file.name = "Na_RINT_CMplot_Manhattan",  
  width = 14,
  height = 6
)

# ===========================
#   QQ plot (qqman)
# ===========================
png(out_qq, width = 1600, height = 1600, res = 200)
qq(gwas$P, main = "QQ Plot — height GWAS")
dev.off()

cat("Manhattan + QQ plots generated in:\n", outdir, "\n")
