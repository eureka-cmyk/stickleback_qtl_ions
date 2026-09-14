library(ggplot2)

# VERY IMPORTANT: avoid X11 fonts
options(bitmapType = "cairo")

ph <- read.table("/gpfs01/home/shxyt1/data/sticklebacks/pheno/traits/pheno_Na.txt",
                 header=FALSE)

colnames(ph) <- c("FID","IID","Na")

out_png <- "/gpfs01/home/shxyt1/data/sticklebacks/pheno/traits/Na_distribution.png"

png(out_png, width=1600, height=1200, res=150)

ggplot(ph, aes(x=Na)) +
  geom_histogram(bins=40, fill="steelblue", color="black") +
  theme_bw() +
  ggtitle("Distribution of Na")

dev.off()

cat("Saved plot to:", out_png, "\n")






# RINT transform
# ph$Na_rint <- qnorm((rank(ph$Na, na.last="keep") - 0.5) / sum(!is.na(ph$Na)))

# 保存 RINT 版表型
write.table(
  ph[, c("FID","IID","Na_rint")],
  "/gpfs01/home/shxyt1/data/sticklebacks/pheno/traits/pheno_Na_RINT.txt",
  row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t"
)


ph_rint <- read.table("/gpfs01/home/shxyt1/data/sticklebacks/pheno/traits/pheno_Na_RINT.txt",
                 header=FALSE)

colnames(ph_rint) <- c("FID","IID","Na")

out_png_RINT <- "/gpfs01/home/shxyt1/data/sticklebacks/pheno/traits/Na_distribution_RINT.png"

png(out_png_RINT, width=1600, height=1200, res=150)

ggplot(ph, aes(x=Na_rint)) +
  geom_histogram(bins=40, fill="steelblue", color="black") +
  theme_bw() +
  ggtitle("Distribution of Na after RINT")

dev.off()

cat("Saved plot to:", out_png_RINT, "\n")



#height_weight_pc5.txt


all_pheno <- "/gpfs01/home/shxyt1/data/sticklebacks/pheno/all_pheno_clean.txt"
pcs_file  <- "/gpfs01/home/shxyt1/data/sticklebacks/pheno/traits/pcs5.txt"


out_file <- "/gpfs01/home/shxyt1/data/sticklebacks/pheno/traits/height_weight_pc5.txt"


ph  <- fread(all_pheno, header = TRUE)
pcs <- fread(pcs_file, header = FALSE)
setnames(pcs, c("FID","IID","PC1","PC2","PC3","PC4","PC5"))

covar <- ph[, .(FID, IID,
                body_length,
                body_weight,
                wet_weight_body_no_head)]

covar <- merge(covar, pcs, by = c("FID", "IID"))

write.table(
  covar,
  file = out_file,
  row.names = FALSE,
  col.names = FALSE,   
  quote = FALSE,
  sep = "\t"
)

cat("Written covariate file (no header) to:\n", out_file, "\n")

