library(data.table)
library(ggplot2)

pca <- fread("/gpfs01/home/shxyt1/data/sticklebacks/VCF/eda/eda_edaRegion.pca.eigenvec",
             col.names = c("FID","IID",
                           paste0("EDA_PC",1:10)))

head(pca)


ggplot(pca, aes(x = EDA_PC1, y = EDA_PC2)) +
  geom_point(size = 2, alpha = 0.7) +
  theme_classic(base_size = 14) +
  labs(
    x = "EDA PC1",
    y = "EDA PC2",
    title = "Local PCA of EDA region"
  )

ggplot(pca, aes(x = EDA_PC2, y = EDA_PC3)) +
  geom_point(size = 2, alpha = 0.7) +
  theme_classic(base_size = 14) +
  labs(
    x = "EDA PC1",
    y = "EDA PC2",
    title = "Local PCA of EDA region"
  )



pca25 <- fread("/gpfs01/home/shxyt1/data/sticklebacks/VCF/eda/eda25kb.pca.eigenvec",
               col.names = c("FID","IID", paste0("EDA50_PC",1:10)))

ggplot(pca25, aes(EDA50_PC1, EDA50_PC2)) +
  geom_point(size = 2, alpha = 0.7) +
  theme_classic(base_size = 14) +
  labs(
    x = "EDA PC1 (±25 kb)",
    y = "EDA PC2 (±25 kb)",
    title = "Local PCA of EDA region (±25 kb)"
  )

pca10 <- fread("/gpfs01/home/shxyt1/data/sticklebacks/VCF/eda/eda10kb.pca.eigenvec",
               col.names = c("FID","IID", paste0("EDA10_PC",1:10)))

ggplot(pca10, aes(EDA10_PC1, EDA10_PC2)) +
  geom_point(size = 2, alpha = 0.7) +
  theme_classic(base_size = 14) +
  labs(
    x = "EDA PC1 (±10 kb)",
    y = "EDA PC2 (±10 kb)",
    title = "Local PCA of EDA region (±10 kb)"
  )


pcag <- fread("/gpfs01/home/shxyt1/data/sticklebacks/VCF/eda/eda_gene.pca.eigenvec",
              col.names = c("FID","IID", paste0("EDA_gene_PC",1:10)))

ggplot(pcag, aes(EDA_gene_PC1, EDA_gene_PC2)) +
  geom_point(size = 2, alpha = 0.7) +
  theme_classic(base_size = 14) +
  labs(
    x = "EDA PC1 (gene body)",
    y = "EDA PC2 (gene body)",
    title = "Local PCA of EDA gene body"
  )



library(data.table)

setwd("/gpfs01/home/shxyt1/data/sticklebacks/VCF/eda")

read_pca2 <- function(file, prefix){
  p <- fread(file)
  setnames(p, c("FID","IID", paste0(prefix, "_PC", 1:(ncol(p)-2))))
  p[, .(IID, PC1 = get(paste0(prefix, "_PC1")), PC2 = get(paste0(prefix, "_PC2")))]
}

assign_class <- function(dt, pc1_col){
  pc1 <- dt[[pc1_col]]
  qs  <- quantile(pc1, probs = c(1/3, 2/3), na.rm = TRUE)
  dt[, c("t1","t2") := .(qs[1], qs[2])]
  dt[, class := fifelse(is.na(get(pc1_col)), NA_integer_,
                        fifelse(get(pc1_col) <= t1, 0L,
                                fifelse(get(pc1_col) <= t2, 1L, 2L)))]
  dt
}


p100 <- read_pca2("eda_edaRegion.pca.eigenvec", "EDA100")
p50  <- read_pca2("eda50kb.pca.eigenvec",      "EDA50")
p25  <- read_pca2("eda25kb.pca.eigenvec",      "EDA25")
p10  <- read_pca2("eda10kb.pca.eigenvec",      "EDA10")
pg   <- read_pca2("eda_gene.pca.eigenvec",     "EDAgene")


setnames(p100, c("PC1","PC2"), c("EDA100_PC1","EDA100_PC2"))
setnames(p50,  c("PC1","PC2"), c("EDA50_PC1","EDA50_PC2"))
setnames(p25,  c("PC1","PC2"), c("EDA25_PC1","EDA25_PC2"))
setnames(p10,  c("PC1","PC2"), c("EDA10_PC1","EDA10_PC2"))
setnames(pg,   c("PC1","PC2"), c("EDAgene_PC1","EDAgene_PC2"))


eda_all <- Reduce(function(x,y) merge(x,y,by="IID",all=TRUE),
                  list(p100,p50,p25,p10,pg))

tmp <- copy(eda_all)
tmp <- assign_class(tmp, "EDA100_PC1"); eda_all[, class_100 := tmp$class]
tmp <- assign_class(tmp, "EDA50_PC1");  eda_all[, class_50  := tmp$class]
tmp <- assign_class(tmp, "EDA25_PC1");  eda_all[, class_25  := tmp$class]
tmp <- assign_class(tmp, "EDA10_PC1");  eda_all[, class_10  := tmp$class]
tmp <- assign_class(tmp, "EDAgene_PC1");eda_all[, class_gene:= tmp$class]


out <- eda_all[, .(
  IID,
  EDA100_PC1, EDA100_PC2, class_100,
  EDA50_PC1,  EDA50_PC2,  class_50,
  EDA25_PC1,  EDA25_PC2,  class_25,
  EDA10_PC1,  EDA10_PC2,  class_10,
  EDAgene_PC1,EDAgene_PC2,class_gene
)]

fwrite(out, "EDA_PCA_all_windows_with_classes.txt", sep="\t")



# =========================
# EDA local PCA classification (1D k-means on PC1)
# + publication-ready colors
# + optional reproducible reassignment for borderline samples
# =========================

library(data.table)
library(ggplot2)

# ---- User settings ----
setwd("/gpfs01/home/shxyt1/data/sticklebacks/VCF/eda")

# Input file: ±25 kb PCA eigenvec (PLINK output)
pca_file <- "eda25kb.pca.eigenvec"

# Output figure names
out_pdf <- "EDA_PCA_25kb_1Dkmeans.pdf"
out_png <- "EDA_PCA_25kb_1Dkmeans.png"

# Random seed for reproducibility
set.seed(123)

# Okabe–Ito colorblind-friendly palette (good for journals)
# Blue, orange, green, etc. Here we use:
# cc (0) = blue, cl (1) = orange, ll (2) = vermillion
pal <- c("0" = "#0072B2",   # blue
         "1" = "#E69F00",   # orange
         "2" = "#D55E00")   # vermillion

# ---- 1) Read PCA eigenvec ----
# PLINK eigenvec format: FID IID PC1 PC2 PC3 ...
pca <- fread(pca_file, header = FALSE)

# Rename columns
# V1=FID, V2=IID, V3=PC1, V4=PC2 ...
setnames(pca, c("FID", "IID", paste0("PC", 1:(ncol(pca) - 2))))

# Keep only PC1 and PC2 for plotting/classification
dt <- pca[, .(IID, PC1 = PC1, PC2 = PC2)]

# ---- 2) 1D k-means on PC1 ----
km <- kmeans(dt[, .(PC1)], centers = 3, nstart = 50)
dt[, cluster := km$cluster]

# Map clusters to ordered classes 0/1/2 by mean(PC1)
cent <- dt[, .(mean_PC1 = mean(PC1)), by = cluster][order(mean_PC1)]
cent[, EDA_class_1D := 0:2]

dt <- merge(dt, cent[, .(cluster, EDA_class_1D)], by = "cluster", all.x = TRUE)
dt[, EDA_class_1D := as.integer(EDA_class_1D)]  # ensure integer 0/1/2

# ---- 3) Identify borderline samples (optional) ----
# We define borderline samples as those whose PC1 is near the class boundaries.
# This helps explain occasional "odd" points without manual editing.
# Boundary positions are midpoints between the ordered class mean(PC1).
means <- cent[, .(EDA_class_1D, mean_PC1)]
setorder(means, EDA_class_1D)

b01 <- mean(means$mean_PC1[means$EDA_class_1D == 0],
            means$mean_PC1[means$EDA_class_1D == 1])
b12 <- mean(means$mean_PC1[means$EDA_class_1D == 1],
            means$mean_PC1[means$EDA_class_1D == 2])

# Define a small window around boundaries to flag uncertainty
# You can adjust "delta" if you want stricter/looser flagging.
delta <- 0.005
dt[, borderline := (abs(PC1 - b01) < delta) | (abs(PC1 - b12) < delta)]

# ---- 4) OPTIONAL: reproducible reassignment for suspicious points ----
# Compute 2D centroids for each class based on current assignment
cent2d <- dt[, .(cPC1 = mean(PC1), cPC2 = mean(PC2)), by = EDA_class_1D]

# Extract centroid coordinates as scalars (avoids merge issues)
c0 <- cent2d[EDA_class_1D == 0, .(cPC1, cPC2)]
c1 <- cent2d[EDA_class_1D == 1, .(cPC1, cPC2)]
c2 <- cent2d[EDA_class_1D == 2, .(cPC1, cPC2)]

# Safety check: ensure all centroids exist
stopifnot(nrow(c0) == 1, nrow(c1) == 1, nrow(c2) == 1)

# Squared distance to a centroid
dist2 <- function(x1, x2, c1, c2) (x1 - c1)^2 + (x2 - c2)^2

# Compute distances to each class centroid
dt[, d0 := dist2(PC1, PC2, c0$cPC1, c0$cPC2)]
dt[, d1 := dist2(PC1, PC2, c1$cPC1, c1$cPC2)]
dt[, d2 := dist2(PC1, PC2, c2$cPC1, c2$cPC2)]

# Automatic reassignment rule (reproducible):
# Only reassign class 1 -> class 0 if:
# - PC1 is clearly left of the 0/1 boundary (b01 - delta)
# - and the point is closer to class 0 centroid than class 1 centroid
dt[, reassigned := FALSE]
dt[EDA_class_1D == 1 & PC1 < (b01 - delta) & d0 < d1,
   `:=`(EDA_class_1D = 0L, reassigned = TRUE)]

# Note: If you do NOT want any reassignment, comment out the line above.
# The recommended scientific default is to keep original k-means result and
# use "borderline" flag instead of reassignment.

# ---- 5) Prepare labels and factor for plotting ----
dt[, class_label := factor(EDA_class_1D,
                           levels = c(0, 1, 2),
                           labels = c("cc (0)", "cl (1)", "ll (2)"))]

# ---- 6) Plot (publication-ready) ----
# Tips for journals:
# - Use small points with alpha
# - Keep clean theme
# - Export as PDF for vector graphics + PNG for quick view

p <- ggplot(dt, aes(x = PC1, y = PC2, color = factor(EDA_class_1D))) +
  geom_point(size = 2.0, alpha = 0.85) +
  scale_color_manual(values = pal,
                     breaks = c("0", "1", "2"),
                     labels = c("cc (0)", "cl (1)", "ll (2)")) +
  theme_classic(base_size = 14) +
  theme(
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 11)
  ) +
  labs(
    title = "Local PCA of EDA region (±25 kb)",
    x = "EDA PC1 (±25 kb)",
    y = "EDA PC2 (±25 kb)",
    color = "EDA class (1D k-means)"
  )

# Save figures
ggsave(out_pdf, plot = p, width = 7.5, height = 5.0, units = "in", dpi = 300)
ggsave(out_png, plot = p, width = 7.5, height = 5.0, units = "in", dpi = 300)

# ---- 7) Export a summary table (optional but recommended) ----
# This table helps you keep a record of final class calls and uncertainty flags.
out_tab <- dt[, .(IID, PC1, PC2,
                  EDA_class_1D,
                  borderline,
                  reassigned)]
fwrite(out_tab, "EDA25kb_PCA_1Dkmeans_class_calls.txt", sep = "\t")

# Print quick diagnostics
cat("Class counts:\n")
print(table(dt$EDA_class_1D))

cat("\nNumber of borderline samples:\n")
print(sum(dt$borderline, na.rm = TRUE))

cat("\nNumber of reassigned samples:\n")
print(sum(dt$reassigned, na.rm = TRUE))

cat("\nBoundaries (PC1 midpoints):\n")
cat("b01 =", b01, " b12 =", b12, "\n")

cc_center <- dt[EDA_class_1D == 0, .(
  cPC1 = mean(PC1),
  cPC2 = mean(PC2)
)]

# blue-side yellow candidates
dt[EDA_class_1D == 1][order(PC1)][1:5, .(IID, PC1, PC2)]

dt[IID == "QTL22188", EDA_class_1D := 0L]
dt[IID == "QTL22188", reassigned := TRUE]

dt[IID == "QTL22188", .(IID, EDA_class_1D, PC1, PC2, reassigned)]

library(ggplot2)

ggplot(
  dt,
  aes(x = PC1, y = PC2, color = factor(EDA_class_1D))
) +
  geom_point(size = 2, alpha = 0.85) +
  scale_color_manual(
    values = c(
      "0" = "#0072B2",  # cc (blue)
      "1" = "#E69F00",  # cl (orange)
      "2" = "#D55E00"   # ll (red)
    ),
    labels = c("cc (0)", "cl (1)", "ll (2)")
  ) +
  theme_classic() +
  labs(
    x = "EDA PC1 (±25 kb)",
    y = "EDA PC2 (±25 kb)",
    color = "EDA class"
  )


# Save figures (recommended for journals)
ggsave("EDA_PCA_25kb_1Dkmeans.pdf", plot = p,
       width = 7.5, height = 5.0, units = "in", dpi = 300)

ggsave("EDA_PCA_25kb_1Dkmeans.png", plot = p,
       width = 7.5, height = 5.0, units = "in", dpi = 300)

# Also print to the RStudio plot window
print(p)


library(data.table)

library(data.table)

run_eda_kmeans_1d <- function(eigenvec_file, window_label) {
  
  # read PCA eigenvec
  pca <- fread(eigenvec_file, header = FALSE)
  setnames(pca, c("FID", "IID", paste0("PC", seq_len(ncol(pca) - 2))))
  
  dt <- pca[, .(IID, PC1, PC2)]
  
  # 1D k-means on PC1
  set.seed(123)
  km <- kmeans(dt[, .(PC1)], centers = 3, nstart = 50)
  dt[, cluster := km$cluster]
  
  # order clusters by mean PC1 -> 0 / 1 / 2
  cent <- dt[, .(mean_PC1 = mean(PC1)), by = cluster][order(mean_PC1)]
  cent[, class := 0:2]
  
  dt <- merge(dt, cent[, .(cluster, class)], by = "cluster", all.x = TRUE)
  
  # ----- SAFE output construction (NO := , NO NSE) -----
  out <- data.table(
    IID = dt$IID
  )
  
  out[[paste0(window_label, "_PC1")]]   <- dt$PC1
  out[[paste0(window_label, "_PC2")]]   <- dt$PC2
  out[[paste0(window_label, "_class")]] <- dt$class
  
  return(out)
}


res_100  <- run_eda_kmeans_1d("eda_edaRegion.pca.eigenvec", "EDA100")
res_50   <- run_eda_kmeans_1d("eda50kb.pca.eigenvec",  "EDA50")
res_25   <- run_eda_kmeans_1d("eda25kb.pca.eigenvec",  "EDA25")
res_10   <- run_eda_kmeans_1d("eda10kb.pca.eigenvec",  "EDA10")
res_gene <- run_eda_kmeans_1d("eda_gene.pca.eigenvec", "EDAgene")

# merge all results by IID
eda_kmeans_all <- Reduce(
  function(x, y) merge(x, y, by = "IID", all = TRUE),
  list(res_100, res_50, res_25, res_10, res_gene)
)

any(duplicated(colnames(eda_kmeans_all)))
colnames(eda_kmeans_all)

fwrite(
  eda_kmeans_all,
  "EDA_kmeans_1D_all_windows.txt",
  sep = "\t"
)


library(data.table)
library(ggplot2)

plot_eda_pca_kmeans <- function(
    eigenvec_file,     # e.g. "eda25kb.pca.eigenvec"
    window_label,      # e.g. "EDA25"
    out_prefix         # e.g. "EDA25"
) {
  # read PCA eigenvec
  pca <- fread(eigenvec_file, header = FALSE)
  setnames(pca, c("FID", "IID", paste0("PC", seq_len(ncol(pca) - 2))))
  dt <- pca[, .(IID, PC1, PC2)]
  
  # 1D k-means on PC1
  set.seed(123)
  km <- kmeans(dt[, .(PC1)], centers = 3, nstart = 50)
  dt[, cluster := km$cluster]
  
  # order clusters by mean PC1 -> 0 / 1 / 2
  cent <- dt[, .(mean_PC1 = mean(PC1)), by = cluster][order(mean_PC1)]
  cent[, class := 0:2]
  dt <- merge(dt, cent[, .(cluster, class)], by = "cluster", all.x = TRUE)
  
  # plot (journal-friendly colors)
  p <- ggplot(
    dt,
    aes(x = PC1, y = PC2, color = factor(class))
  ) +
    geom_point(size = 2, alpha = 0.85) +
    scale_color_manual(
      values = c(
        "0" = "#0072B2",  # cc
        "1" = "#E69F00",  # cl
        "2" = "#D55E00"   # ll
      ),
      labels = c("cc (0)", "cl (1)", "ll (2)")
    ) +
    theme_classic(base_size = 14) +
    labs(
      title = paste0("Local PCA of EDA region (", window_label, ")"),
      x = paste0("EDA PC1 (", window_label, ")"),
      y = paste0("EDA PC2 (", window_label, ")"),
      color = "EDA class (1D k-means)"
    )
  
  # save figures
  ggsave(paste0(out_prefix, "_PCA_kmeans.pdf"),
         plot = p, width = 10, height = 5, units = "in", dpi = 300)
  ggsave(paste0(out_prefix, "_PCA_kmeans.png"),
         plot = p, width = 10, height = 5, units = "in", dpi = 300)
  
  return(invisible(dt))
}

plot_eda_pca_kmeans("eda_edaRegion.pca.eigenvec", "±100 kb", "EDA100")
plot_eda_pca_kmeans("eda50kb.pca.eigenvec",       "±50 kb",  "EDA50")
plot_eda_pca_kmeans("eda25kb.pca.eigenvec",       "±25 kb",  "EDA25")
plot_eda_pca_kmeans("eda10kb.pca.eigenvec",       "±10 kb",  "EDA10")
plot_eda_pca_kmeans("eda_gene.pca.eigenvec",      "gene",    "EDAgene")

