# Set working directory (change to your actual path)
# setwd("/gpfs01/home/shxyt1/data/sticklebacks/VCF/eda/")

"/gpfs01/home/shxyt1/data/sticklebacks/VCF/eda/"

# Load the sample missingness file
# Header: FID, IID, MISS_PHENO, N_MISS, N_GENO, F_MISS
#sample
imiss <- read.table("/gpfs01/home/shxyt1/data/sticklebacks/VCF/eda/sample_missing_10kb.imiss", header = TRUE, stringsAsFactors = FALSE)

#snp
lmiss <- read.table("/gpfs01/home/shxyt1/data/sticklebacks/VCF/eda/sample_missing_10kb.lmiss", header = TRUE, stringsAsFactors = FALSE)

# 1. Calculate basic descriptive statistics for missing rate
# This helps you understand the overall data quality
missing_summary <- summary(imiss$F_MISS)
print("Summary of Missing Rate (F_MISS):")
print(missing_summary)

# 2. Identify samples with missing rate > 10%
# You can adjust '0.1' to your preferred threshold
high_missing_samples <- imiss[imiss$F_MISS > 0.1, ]
print("Samples with missing rate > 10%:")
print(high_missing_samples)

# 3. Create a histogram to visualize the distribution
library(ggplot2)
ggplot(lmiss, aes(x = F_MISS)) +
  geom_histogram(binwidth = 0.01, fill = "steelblue", color = "white") +
  geom_vline(xintercept = 0.1, linetype = "dashed", color = "red") + # Threshold line
  labs(title = "Sample Missingness Distribution",
       subtitle = "Red line indicates 10% missing rate threshold",
       x = "Missing Rate (F_MISS)",
       y = "Number of Samples") +
  theme_minimal()







# Load the sample missingness file
# Header: FID, IID, MISS_PHENO, N_MISS, N_GENO, F_MISS
lmiss <- read.table("/gpfs01/home/shxyt1/data/sticklebacks/VCF/eda/sample_missing_10kb.lmiss", header = TRUE, stringsAsFactors = FALSE)

# 1. Calculate basic descriptive statistics for missing rate
# This helps you understand the overall data quality
missing_summary <- summary(lmiss$F_MISS)
print("Summary of Missing Rate (F_MISS):")
print(missing_summary)

# 2. Identify samples with missing rate > 10%
# You can adjust '0.1' to your preferred threshold
high_missing_samples <- lmiss[lmiss$F_MISS > 0.1, ]
print("Samples with missing rate > 10%:")
print(high_missing_samples)

# 3. Create a histogram to visualize the distribution
library(ggplot2)
ggplot(lmiss, aes(x = F_MISS)) +
  geom_histogram(binwidth = 0.01, fill = "steelblue", color = "white") +
  geom_vline(xintercept = 0.1, linetype = "dashed", color = "red") + # Threshold line
  labs(title = "SNP Missingness Distribution",
       subtitle = "Red line indicates 10% missing rate threshold",
       x = "Missing Rate (F_MISS)",
       y = "Number of Samples") +
  theme_minimal()









library(data.table)
library(ggplot2)

# =========================
# File paths
# =========================
base_path <- "/gpfs01/home/shxyt1/data/sticklebacks/VCF/eda"

imiss_file <- file.path(base_path, "eda10kb_geno0.01.filt.named.imiss")
lmiss_file <- file.path(base_path, "eda10kb_geno0.01.filt.named.lmiss")
raw_file   <- file.path(base_path, "eda10kb_geno0.01.filt.named.raw")
pca_file   <- file.path(base_path, "eda10kb_geno0.01.filt.named.eigenvec")
mds_file   <- file.path(base_path, "eda10kb_geno0.01.filt.named.mds")
bim_file   <- file.path(base_path, "eda10kb_geno0.01.filt.named.bim")

out_dir <- file.path(base_path, "missingness_final_plots")
dir.create(out_dir, showWarnings = FALSE)

# =========================
# Read sample-level data
# =========================
imiss <- fread(imiss_file)

pca <- fread(pca_file, header = FALSE)
setnames(pca, c("FID","IID",paste0("PC",1:(ncol(pca)-2))))
pca <- pca[, .(FID, IID, PC1)]

mds <- fread(mds_file)
mds <- mds[, .(FID, IID, MDS1 = C1)]

sample_dt <- merge(imiss[, .(FID, IID, F_MISS)], pca, by=c("FID","IID"))
sample_dt <- merge(sample_dt, mds, by=c("FID","IID"))

# =========================
# Plot 1: sample vs PC1
# =========================
p1 <- ggplot(sample_dt, aes(PC1, F_MISS)) +
  geom_point(alpha=0.7) +
  geom_smooth(method="lm") +
  theme_bw() +
  labs(title="Sample missing vs PC1")
p1
ggsave(file.path(out_dir,"sample_vs_PC1-geno0.01-mind0.1.png"), p1, width=6, height=5)

# =========================
# Plot 2: sample vs MDS1
# =========================
p2 <- ggplot(sample_dt, aes(MDS1, F_MISS)) +
  geom_point(alpha=0.7) +
  geom_smooth(method="lm") +
  theme_bw() +
  labs(title="Sample missing vs MDS1")
p2
ggsave(file.path(out_dir,"sample_vs_MDS1-geno0.01-mind0.1.png"), p2, width=6, height=5)

# =========================
# Read genotype matrix
# =========================
raw_dt <- fread(raw_file)

meta_cols <- c("FID","IID","PAT","MAT","SEX","PHENOTYPE")
geno_cols <- setdiff(colnames(raw_dt), meta_cols)

geno <- raw_dt[, ..geno_cols]

# =========================
# Add PC1 & MDS1
# =========================
geno_merge <- cbind(geno, PC1 = sample_dt$PC1, MDS1 = sample_dt$MDS1)

# =========================
# Correlation function
# =========================
calc_cor <- function(g, x){
  ok <- !is.na(g) & !is.na(x)
  if(sum(ok) < 3) return(NA)
  if(sd(g[ok]) == 0) return(NA)
  cor(g[ok], x[ok])
}

# =========================
# SNP correlations
# =========================
pc1_cor  <- sapply(geno_cols, function(s) calc_cor(geno_merge[[s]], geno_merge$PC1))
mds1_cor <- sapply(geno_cols, function(s) calc_cor(geno_merge[[s]], geno_merge$MDS1))

snp_dt <- data.table(
  SNP = geno_cols,
  cor_PC1 = pc1_cor,
  cor_MDS1 = mds1_cor
)

# =========================
# Add missingness (直接按顺序)
# =========================
lmiss <- fread(lmiss_file)

snp_dt[, F_MISS := lmiss$F_MISS]

# =========================
# Plot 3: SNP vs PC1
# =========================
p3 <- ggplot(snp_dt, aes(cor_PC1, F_MISS)) +
  geom_point(alpha=0.7) +
  geom_smooth(method="lm") +
  theme_bw() +
  labs(title="SNP missing vs PC1 correlation")

ggsave(file.path(out_dir,"SNP_vs_PC1-geno0.01-mind0.1.png"), p3, width=6, height=5)

# =========================
# Plot 4: SNP vs MDS1
# =========================
p4 <- ggplot(snp_dt, aes(cor_MDS1, F_MISS)) +
  geom_point(alpha=0.7) +
  geom_smooth(method="lm") +
  theme_bw() +
  labs(title="SNP missing vs MDS1 correlation")

ggsave(file.path(out_dir,"SNP_vs_MDS1-geno0.01-mind0.1.png"), p4, width=6, height=5)

cor(sample_dt$PC1, sample_dt$MDS1)
# =========================
# Print correlations
# =========================
cat("Sample vs PC1:", cor(sample_dt$F_MISS, sample_dt$PC1), "\n")
cat("Sample vs MDS1:", cor(sample_dt$F_MISS, sample_dt$MDS1), "\n")

cat("SNP vs PC1:", cor(snp_dt$F_MISS, snp_dt$cor_PC1, use="complete.obs"), "\n")
cat("SNP vs MDS1:", cor(snp_dt$F_MISS, snp_dt$cor_MDS1, use="complete.obs"), "\n")
