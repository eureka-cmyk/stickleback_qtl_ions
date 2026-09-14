library(data.table)
library(ggplot2)
library(scales)

options(bitmapType = "cairo")

# =========================================
# 1. Input / output
# =========================================
infile  <- "/gpfs01/home/shxyt1/data/sticklebacks/pheno/all_pheno_clean.txt"
out_png <- "/gpfs01/home/shxyt1/data/sticklebacks/pheno/elements_boxplot_piecewise_axis_sorted_by_lowerwhisker.png"
out_pdf <- "/gpfs01/home/shxyt1/data/sticklebacks/pheno/elements_boxplot_piecewise_axis_sorted_by_lowerwhisker.pdf"

# =========================================
# 2. Filtering thresholds
# =========================================
missing_threshold <- 0.5

# =========================================
# 3. Read data
# =========================================
dt <- fread(infile, sep = "\t", header = TRUE, data.table = FALSE)

# =========================================
# 4. Element traits only
# =========================================
element_traits <- c(
  "Na","Mg","P","S","K","Ca","Al","V","Mn","Fe",
  "Co","Cu","Zn","As","Se","Rb","Sr","Mo","Ag","Cd","Cs","Ba","Tl","Pb","U"
)

element_traits <- intersect(element_traits, colnames(dt))

# Convert to numeric
for (v in element_traits) {
  dt[[v]] <- suppressWarnings(as.numeric(dt[[v]]))
}

# =========================================
# 5. Remove traits with too many missing values
# =========================================
trait_qc <- rbindlist(lapply(element_traits, function(v) {
  x <- dt[[v]]
  x_non_na <- x[!is.na(x)]
  
  data.table(
    trait = v,
    missing_fraction = mean(is.na(x)),
    negative_fraction = if (length(x_non_na) > 0) mean(x_non_na < 0) else NA_real_,
    non_na_n = sum(!is.na(x)),
    unique_non_na_n = length(unique(x_non_na))
  )
}))

negative_threshold <- 0.5

element_traits_keep <- trait_qc[
  missing_fraction <= missing_threshold &
    negative_fraction <= negative_threshold &
    non_na_n > 1 &
    unique_non_na_n > 1,
  trait
]

print(trait_qc[order(-negative_fraction)])


cat("Traits retained after missing-value filtering:\n")
print(element_traits_keep)

# =========================================
# 6. Reshape to long format
# =========================================
long_dt <- melt(
  as.data.table(dt[, element_traits_keep, drop = FALSE]),
  measure.vars = element_traits_keep,
  variable.name = "trait",
  value.name = "value"
)

# Keep only finite non-missing values and remove negatives
long_dt <- long_dt[!is.na(value) & is.finite(value) & value >= 0]

# =========================================
# 7. Remove traits with too few values after filtering
# =========================================
trait_qc_after <- long_dt[, .(
  n_after_filter = .N,
  unique_after_filter = uniqueN(value)
), by = trait]

traits_final <- trait_qc_after[
  n_after_filter > 1 & unique_after_filter > 1,
  trait
]

long_dt <- long_dt[trait %in% traits_final]

cat("Traits retained after removing negatives:\n")
print(sort(unique(long_dt$trait)))

# =========================================
# 8. Function to calculate lower whisker
#    lower whisker = smallest value >= Q1 - 1.5*IQR
# =========================================
# =========================================
# 9. Order traits by the bottom of the box (Q1)
# =========================================
trait_order_dt <- long_dt[, .(
  q1 = quantile(value, 0.25, na.rm = TRUE, type = 7)
), by = trait]

trait_order_dt <- trait_order_dt[order(q1, trait)]

cat("Trait order by Q1 (bottom of box):\n")
print(trait_order_dt)

trait_order <- trait_order_dt$trait
long_dt[, trait := factor(trait, levels = trait_order)]

# =========================================
# 10. Low-saturation color palette
# =========================================
trait_colors <- c(
  "Na" = "#9FC6CF",
  "Mg" = "#8FB7A7",
  "P"  = "#8C96B8",
  "S"  = "#D9B7A8",
  "K"  = "#A8AFCA",
  "Ca" = "#B7D4CC",
  "Al" = "#8A9CC2",
  "V"  = "#D09A9A",
  "Mn" = "#A28AB4",
  "Fe" = "#8EB8B2",
  "Co" = "#C48E98",
  "Cu" = "#C090AA",
  "Zn" = "#A9A9A9",
  "As" = "#8E8E8E",
  "Se" = "#D2B089",
  "Rb" = "#B7C7B0",
  "Sr" = "#D9C37A",
  "Mo" = "#B0B0B0",
  "Ag" = "#C79A95",
  "Cd" = "#A4BED7",
  "Cs" = "#7F98AD",
  "Ba" = "#B3A26E",
  "Tl" = "#A9968B",
  "Pb" = "#A58AB6",
  "U"  = "#8FB68D"
)

trait_colors <- trait_colors[trait_order]

# =========================================
# 11. Custom piecewise transformation
#     0 to 1  -> linear, length = 1
#     1 to 100 -> same visual length = 1
#     every 100-fold increase above 1 gets equal spacing
# =========================================
piecewise_trans <- trans_new(
  name = "piecewise_0_1_vs_1_100",
  transform = function(x) {
    ifelse(
      x <= 1,
      x,
      1 + log10(x) / 1.5
      
    )
  },
  inverse = function(x) {
    ifelse(
      x <= 1,
      x,
      10^((x - 1) * 2)
    )
  }
)

# =========================================
# 12. Axis breaks
# =========================================
y_breaks <- c(0, 0.25, 0.5, 0.75, 1, 10, 100, 1000, 10000,50000)
y_labels <- c("0", "0.25", "0.5", "0.75", "1", "10", "100", "1k", "10k", "50k")

# =========================================
# 13. Plot
# =========================================
p <- ggplot(long_dt, aes(x = trait, y = value, color = trait)) +
  geom_boxplot(
    width = 0.58,
    fill = NA,
    linewidth = 0.65,
    coef = 4,
    outlier.shape = 16,
    outlier.size = 1.8,
    outlier.alpha = 0.85
  ) +
  scale_color_manual(values = trait_colors, drop = FALSE) +
  scale_y_continuous(
    trans = piecewise_trans,
    breaks = y_breaks,
    labels = y_labels,
    limits = c(0, NA),
    expand = c(0, 0)
  ) +
  labs(
    x = NULL,
    y = "Concentration"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,    
      vjust = 1,    
      size = 11,
      color = "black"
    ),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title.y = element_text(size = 12),
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.line.y = element_line(color = "black", linewidth = 0.5),
    legend.position = "none",
    plot.margin = margin(10, 10, 10, 10)
  )

print(p)

# =========================================
# 14. Save
# =========================================
ggsave(out_png, p, width = 11, height = 6, dpi = 500)
ggsave(out_pdf, p, width = 16, height = 6)

cat("Saved to:\n")
cat(out_png, "\n")
cat(out_pdf, "\n")

