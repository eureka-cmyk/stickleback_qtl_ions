library(data.table)
options(bitmapType = "cairo")

het <- fread("/gpfs01/home/shxyt1/data/sticklebacks/VCF/eda/eda25kb.het")

het[, HET_RATE := 1 - `O(HOM)` / `N(NM)`]

head(het)

fwrite(
  het[, .(IID, HET_RATE)],
  "/gpfs01/home/shxyt1/data/sticklebacks/VCF/eda/eda25kb.heterozygosity.txt",
  sep = "\t"
)


library(data.table)
library(ggplot2)

# =========================
# 0) Paths (edit if needed)
# =========================
# eigenvec_path <- "/gpfs01/home/shxyt1/data/sticklebacks/VCF/eda/eda25kb.pca.eigenvec"
eigenvec_path <- "/gpfs01/home/shxyt1/data/sticklebacks/VCF/eda/eda10kb.pca.eigenvec"

# het_path      <- "/gpfs01/home/shxyt1/data/sticklebacks/VCF/eda/eda25kb.het"
het_path      <- "/gpfs01/home/shxyt1/data/sticklebacks/VCF/eda/eda10kb.het"

# out_het_png   <- "/gpfs01/home/shxyt1/data/sticklebacks/VCF/eda/eda25kb_PC1_vs_Heterozygosity_kmeans.png"
out_het_png   <- "/gpfs01/home/shxyt1/data/sticklebacks/VCF/eda/eda10kb_PC1_vs_Heterozygosity.png"

# out_f_png     <- "/gpfs01/home/shxyt1/data/sticklebacks/VCF/eda/eda25kb_PC1_vs_F_kmeans.png"
out_f_png     <- "/gpfs01/home/shxyt1/data/sticklebacks/VCF/eda/eda10kb_PC1_vs_F.png"

#window_label  <- "12,783,500-12,800,300"

# =========================
# 1) Read PCA (PC1)
# =========================
pca <- fread(eigenvec_path, header = FALSE)
# Standard PLINK: FID IID PC1 PC2 ...
setnames(pca, c("FID", "IID", paste0("PC", 1:(ncol(pca)-2))))
pca <- pca[, .(IID, PC1)]

# =========================
# 2) 1D k-means on PC1 (k=3)
#    Map clusters by mean PC1 -> cc(0), cl(1), ll(2)
# =========================
set.seed(1)

km <- kmeans(pca$PC1, centers = 3, nstart = 50)
pca[, km_raw := km$cluster]

# cluster means (3 rows)
centers_mean <- pca[, .(meanPC1 = mean(PC1)), by = km_raw][order(meanPC1)]

# robust mapping: rank clusters by mean PC1 -> 0,1,2
centers_mean[, class_id := seq_len(.N) - 1L]

# keep only mapping columns (ensure uniqueness)
map <- unique(centers_mean[, .(km_raw, class_id)])

# join back (safe one-to-one)
pca <- merge(pca, map, by = "km_raw", all.x = TRUE, allow.cartesian = FALSE)


pca[, class_lab := factor(
  class_id,
  levels = c(0, 1, 2),
  labels = c("0", "1", "2")
)]

# Colors to match your template style (blue / gold / orange)
cols <- c(
  "0" = "#1f77b4",
  "1" = "#E6AB02",
  "2" = "#D95F02"
)

# =========================
# 3) Read PLINK .het and compute true heterozygosity
# =========================
het <- fread(het_path)

# True heterozygosity rate: 1 - O(HOM)/N(NM)
het[, HET_RATE := 1 - `O(HOM)` / `N(NM)`]

het <- het[, .(IID, F, HET_RATE)]

# =========================
# 4) Merge
# =========================
dat <- merge(pca, het, by = "IID")
cat("Samples after merge:", nrow(dat), "\n")

# =========================
# 5) Plot template function
# =========================
make_template_plot <- function(df, yvar, ylab, title_txt) {
  ggplot(df, aes(x = PC1, y = .data[[yvar]], color = class_lab)) +
    geom_point(size = 2, alpha = 0.85) +
    scale_color_manual(values = cols, name = "EDA class (1D k-means)") +
    theme_classic(base_size = 14)+
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "right",
      
      panel.grid = element_blank(),     # remove grid
      axis.line = element_line(color="black"),
      axis.ticks = element_line(color="black")
    )+
    labs(
      title = title_txt,
      x = paste0("EDA PC1"), 
      y = ylab
    )
}

# =========================
# 6) Plot 1: PC1 vs Heterozygosity
# =========================
p_het <- make_template_plot(
  dat,
  yvar = "HET_RATE",
  ylab = "Heterozygosity",
  title_txt = "Local PCA of EDA region (±10kb)"
)

print(p_het)
ggsave(out_het_png, plot = p_het, width = 9, height = 5, dpi = 300)
cat("Saved:", out_het_png, "\n")

# =========================
# 7) Plot 2: PC1 vs F coefficient
# =========================
p_f <- make_template_plot(
  dat,
  yvar = "F",
  ylab = "Inbreeding coefficient (F)",
  title_txt = "Local PCA of EDA region (±25 kb)"
)

print(p_f)
ggsave(out_f_png, plot = p_f, width = 7.5, height = 5, dpi = 300)
cat("Saved:", out_f_png, "\n")



library(data.table)

# ---- master (summary) table path ----
sum_path <- "/gpfs01/home/shxyt1/data/sticklebacks/VCF/eda/EDA_PCA_all_windows_with_classes.txt"
out_sum  <- "/gpfs01/home/shxyt1/data/sticklebacks/VCF/eda/EDA_PCA_all_windows_with_classes.txt"  # overwrite

# 1) read master table
sumdt <- fread(sum_path)

# 2) build the add-on table from your existing `dat`
#    - keep numeric class_25 consistent with your plot: class_id (0/1/2)
#    - also keep label if you want (optional)
# add <- dat[, .(
#   IID,
#   HET_RATE_sub = HET_RATE,
#   class_sub    = class_id,
#   F_sub        = F
# )]

add <- dat[, .(
  IID,
  HET_RATE_gene = HET_RATE,   # 修改这里
  class_gene    = class_id,   # 修改这里
  F_gene        = F           # 修改这里
)]
# 3) merge into master by IID
sumdt <- merge(sumdt, add, by = "IID", all.x = TRUE)

# 4) insert the new columns right after class_25
cols <- names(sumdt)
pos  <- match("class_25", cols)

# remove existing versions to avoid duplicates if you rerun
base_cols <- cols[!cols %in% c("HET_RATE","hetero_class","F","F_class")]

new_order <- append(base_cols, c("HET_RATE","hetero_class","F","F_class"), after = pos)
setcolorder(sumdt, new_order)

# 5) write out
fwrite(sumdt, out_sum, sep = "\t", quote = FALSE)

cat("Updated summary table written to:\n", out_sum, "\n")






library(data.table)

# Function to process all specific classifications requested
process_eda_final_logic <- function(pca_file, het_file, suffix) {
  
  # 1. Read PCA (PC1, PC2)
  pca <- fread(pca_file, header = FALSE)
  setnames(pca, c("FID", "IID", paste0("PC", 1:(ncol(pca)-2))))
  dt_pca <- pca[, .(IID, PC1, PC2)]
  
  # 2. Read Heterozygosity (.het)
  het <- fread(het_file)
  het[, HET_RATE := 1 - `O(HOM)` / `N(NM)`]
  dt_het <- het[, .(IID, F, HET_RATE)]
  
  # 3. Merge Data
  combined <- merge(dt_pca, dt_het, by = "IID")
  
  # --- [1] 1D K-means on PC1 only ---
  set.seed(1)
  km_pc1 <- kmeans(combined$PC1, centers = 3, nstart = 50)
  combined[, cluster_pc1 := km_pc1$cluster]
  ord_pc1 <- combined[, .(m = mean(PC1)), by = cluster_pc1][order(m)]
  ord_pc1[, PCA_1D_class := 0:2]
  combined <- merge(combined, ord_pc1[, .(cluster_pc1, PCA_1D_class)], by = "cluster_pc1")
  
  # --- [2] 2D K-means on PC1 AND PC2 (Spatial PCA clustering) ---
  scaled_pca <- scale(combined[, .(PC1, PC2)])
  set.seed(1)
  km_pca2d <- kmeans(scaled_pca, centers = 3, nstart = 50)
  combined[, cluster_pca2d := km_pca2d$cluster]
  ord_pca2d <- combined[, .(m = mean(PC1)), by = cluster_pca2d][order(m)]
  ord_pca2d[, PCA_2D_class := 0:2]
  combined <- merge(combined, ord_pca2d[, .(cluster_pca2d, PCA_2D_class)], by = "cluster_pca2d")
  
  # --- [3] 1D K-means on HET_RATE only ---
  set.seed(1)
  km_het <- kmeans(combined$HET_RATE, centers = 3, nstart = 50)
  combined[, cluster_het := km_het$cluster]
  ord_het <- combined[, .(m = mean(HET_RATE)), by = cluster_het][order(m)]
  ord_het[, HET_1D_class := 0:2]
  combined <- merge(combined, ord_het[, .(cluster_het, HET_1D_class)], by = "cluster_het")
  
  # --- [4] Combined 2D K-means on PC1 AND HET_RATE ---
  scaled_pc1_het <- scale(combined[, .(PC1, HET_RATE)])
  set.seed(1)
  km_comb <- kmeans(scaled_pc1_het, centers = 3, nstart = 50)
  combined[, cluster_comb := km_comb$cluster]
  ord_comb <- combined[, .(m = mean(PC1)), by = cluster_comb][order(m)]
  ord_comb[, Combined_PC1_Het_class := 0:2]
  combined <- merge(combined, ord_comb[, .(cluster_comb, Combined_PC1_Het_class)], by = "cluster_comb")
  
  # 4. Final selection and Renaming columns with suffix
  res <- combined[, .(IID, PC1, PC2, HET_RATE, F, 
                      PCA_1D_class, PCA_2D_class, HET_1D_class, Combined_PC1_Het_class)]
  
  colnames(res)[-1] <- paste0(colnames(res)[-1], "_", suffix)
  
  return(res)
}

# --- Execution for the 3 regions ---

# 1. 25kb region
data_25kb <- process_eda_final_logic(
  "/gpfs01/home/shxyt1/data/sticklebacks/VCF/eda/eda25kb.pca.eigenvec",
  "/gpfs01/home/shxyt1/data/sticklebacks/VCF/eda/eda25kb.het",
  "25kb"
)

# 2. 10kb region
data_10kb <- process_eda_final_logic(
  "/gpfs01/home/shxyt1/data/sticklebacks/VCF/eda/eda10kb.pca.eigenvec",
  "/gpfs01/home/shxyt1/data/sticklebacks/VCF/eda/eda10kb.het",
  "10kb"
)

# 3. 16.8kb subset (12.78-12.80 Mb)
data_16kb <- process_eda_final_logic(
  "/gpfs01/home/shxyt1/data/sticklebacks/VCF/eda/eda_subset_16kb.pca.eigenvec",
  "/gpfs01/home/shxyt1/data/sticklebacks/VCF/eda/eda_subset_16kb_het.het",
  "16kb"
)

# --- Final Merge ---
final_eda_summary <- Reduce(function(x, y) merge(x, y, by = "IID", all = TRUE), 
                            list(data_25kb, data_10kb, data_16kb))

# --- Export Result ---
output_file <- "/gpfs01/home/shxyt1/data/sticklebacks/VCF/eda/EDA_Final_Combined_Table.txt"
fwrite(final_eda_summary, output_file, sep = "\t")

# Show in RStudio
View(final_eda_summary)






# 1. Load the trait file
# Since pcs5.txt has no header, we manually name the first two columns as FID and IID
# The remaining columns will be V3, V4, etc.
pheno <- fread("/gpfs01/home/shxyt1/data/sticklebacks/pheno/traits/pcs5.txt", header = FALSE)
setnames(pheno, c("V1", "V2"), c("FID", "IID"))

# 2. Load the EDA summary table
eda_summary <- fread("/gpfs01/home/shxyt1/data/sticklebacks/VCF/eda/EDA_Final_Combined_Table.txt")

# 3. Merge by IID
# We take only the IID and the specific classification column from eda_summary
final_pheno <- merge(pheno, eda_summary[, .(IID, Combined_PC1_Het_class_10kb)], by = "IID", all.x = TRUE)

# 4. Optional: If you want to keep the original PLINK-like order (FID then IID)
# Move FID back to the first column (merge usually puts the 'by' column first)
setcolorder(final_pheno, c("FID", "IID"))

# 5. Export without title (Header)
fwrite(final_pheno, 
       "/gpfs01/home/shxyt1/data/sticklebacks/pheno/traits/pc5_eda_10kb_hetro_class.txt", 
       sep = "\t", 
       col.names = FALSE)

# Check the output in console
head(final_pheno)



#MDS
library(ggplot2)

mds_df <- read.table("/gpfs01/home/shxyt1/data/sticklebacks/VCF/eda/eda10kb_mds.mds", header = TRUE)

# 2. Generate the plot with all points in grey and no categorization
ggplot(mds_df, aes(x = C1, y = C2)) +
  # Set all points to grey, adjust size and transparency
  geom_point(size = 2, alpha = 0.6, color = "black") +
  # Apply classic theme with base font size 14
  theme_classic(base_size = 14) +
  # Fine-tune theme elements to match your template
  theme(
    plot.title = element_text(hjust = 0.5),    # Center the title
    legend.position = "none",                  # Remove legend as there is no classification
    panel.grid = element_blank(),              # Remove all grid lines
    axis.line = element_line(color = "black"),    # Solid black axis lines
    axis.ticks = element_line(color = "black")    # Solid black axis ticks
  ) +
  # Set axis labels and plot title
  labs(
    title = "MDS Plot of EDA Region(±10kb)",
    x = "EDA MDS1", 
    y = "EDA MDS2"
  )
