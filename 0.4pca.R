library(data.table)
library(ggplot2)
library(ggrepel)

options(bitmapType = "cairo")

# =========================================
# 1. Input / output
# =========================================
infile  <- "/gpfs01/home/shxyt1/data/sticklebacks/pheno/all_pheno_clean.txt"
out_png <- "/gpfs01/home/shxyt1/data/sticklebacks/pheno/PCA_variable_contributions_elements.png"
out_pdf <- "/gpfs01/home/shxyt1/data/sticklebacks/pheno/PCA_variable_contributions_elements.pdf"

# =========================================
# 2. Settings
# =========================================
missing_threshold <- 0.5

element_traits <- c(
  "Na","Mg","P","S","K","Ca","Al","V","Mn","Fe",
  "Co","Cu","Zn","As","Se","Rb","Sr","Mo","Ag","Cd","Cs","Ba","Tl","Pb","U"
)
# =========================================
# 3. Read data
# =========================================
dt <- fread(infile, sep = "\t", header = TRUE, data.table = FALSE)

element_traits <- intersect(element_traits, colnames(dt))

for (v in element_traits) {
  dt[[v]] <- suppressWarnings(as.numeric(dt[[v]]))
}

# =========================================
# 4. Filter elements
# =========================================
trait_qc <- rbindlist(lapply(element_traits, function(v) {
  x <- dt[[v]]
  data.table(
    trait = v,
    missing_fraction = mean(is.na(x)),
    non_na_n = sum(!is.na(x)),
    unique_non_na_n = length(unique(x[!is.na(x)]))
  )
}))

element_traits_keep <- trait_qc[
  missing_fraction <= missing_threshold &
    non_na_n > 1 &
    unique_non_na_n > 1,
  trait
]

# Build analysis matrix and remove negative values
pca_dt <- as.data.table(dt[, element_traits_keep, drop = FALSE])

for (v in names(pca_dt)) {
  x <- pca_dt[[v]]
  x[is.na(x) | !is.finite(x)] <- NA_real_
  pca_dt[[v]] <- x
}

# Remove variables with too few values after filtering
keep2 <- names(pca_dt)[sapply(names(pca_dt), function(v) {
  x <- pca_dt[[v]]
  sum(!is.na(x)) > 3 && length(unique(x[!is.na(x)])) > 1
})]

pca_dt <- pca_dt[, ..keep2]

if (ncol(pca_dt) < 2) {
  stop("Fewer than two valid element variables remain for PCA.")
}

# Remove rows with any missing values
pca_complete <- na.omit(pca_dt)

if (nrow(pca_complete) < 5) {
  stop("Too few complete samples remain after filtering for PCA.")
}

cat("Variables used in PCA:\n")
print(colnames(pca_complete))
cat("Number of complete samples used in PCA:\n")
print(nrow(pca_complete))

# =========================================
# 5. Run PCA
# =========================================
pca_res <- prcomp(pca_complete, center = TRUE, scale. = TRUE)

# =========================================
# 6. Build variable coordinates
#    Correlation circle coordinates:
#    loadings * sdev
# =========================================
var_coord <- as.data.table(pca_res$rotation[, 1:2] %*% diag(pca_res$sdev[1:2]))
var_coord[, variable := rownames(pca_res$rotation)]

colnames(var_coord)[1:2] <- c("Dim1", "Dim2")

# Contribution of variables to PC1 + PC2
eig <- pca_res$sdev^2
var_coord[, contrib_PC1 := 100 * (Dim1^2) / sum(Dim1^2)]
var_coord[, contrib_PC2 := 100 * (Dim2^2) / sum(Dim2^2)]
var_coord[, contrib := contrib_PC1 + contrib_PC2]

# Percent variance explained
pc1_pct <- round(100 * eig[1] / sum(eig), 1)
pc2_pct <- round(100 * eig[2] / sum(eig), 1)

# =========================================
# 7. Unit circle
# =========================================
theta <- seq(0, 2 * pi, length.out = 500)
circle_dt <- data.table(
  x = cos(theta),
  y = sin(theta)
)

# =========================================
# 8. Plot
# =========================================
p <- ggplot() +
  geom_path(
    data = circle_dt,
    aes(x = x, y = y),
    color = "grey70",
    linewidth = 0.6
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_segment(
    data = var_coord,
    aes(
      x = 0, y = 0,
      xend = Dim1, yend = Dim2,
      color = contrib
    ),
    arrow = arrow(length = unit(0.18, "cm")),
    linewidth = 0.7
  ) +
  geom_text_repel(
    data = var_coord,
    aes(
      x = Dim1,
      y = Dim2,
      label = variable,
      color = contrib
    ),
    size = 4,
    box.padding = 0.2,
    point.padding = 0.15,
    segment.color = NA,
    show.legend = FALSE
  ) +
  scale_color_gradientn(
    colours = c("#2166AC", "gold", "#B2182B"),
    name = "contrib"
  ) +
  coord_equal(xlim = c(-1.05, 1.05), ylim = c(-1.05, 1.05), expand = FALSE) +
  labs(
    title = "Variable Contributions",
    x = paste0("Dim1 (", pc1_pct, "%)"),
    y = paste0("Dim2 (", pc2_pct, "%)")
  ) +
  # theme_minimal(base_size = 13) +
  theme_classic() +
  theme(
    panel.grid = element_blank(),   # no grid
    axis.line = element_blank(),    # ← remove axis lines
    axis.ticks = element_blank(),   # ← remove ticks
    plot.title = element_text(size = 14),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10)
  )

print(p)

# =========================================
# 9. Save
# =========================================
ggsave(out_png, p, width = 7, height = 6, dpi = 500, bg = "white")
ggsave(out_pdf, p, width = 7, height = 6, bg = "white")

cat("Saved to:\n")
cat(out_png, "\n")
cat(out_pdf, "\n")


range_dt <- data.table(
  element = element_traits,
  min = sapply(element_traits, function(v) min(dt[[v]], na.rm = TRUE)),
  max = sapply(element_traits, function(v) max(dt[[v]], na.rm = TRUE))
)

