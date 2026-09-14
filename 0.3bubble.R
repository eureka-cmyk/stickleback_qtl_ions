library(data.table)
library(ggplot2)

options(bitmapType = "cairo")

# =========================================
# 1. Input / output
# =========================================
infile  <- "/gpfs01/home/shxyt1/data/sticklebacks/pheno/all_pheno_clean.txt"
out_png <- "/gpfs01/home/shxyt1/data/sticklebacks/pheno/element_body_bubble_heatmap_sorted_by_height_rho.png"
out_pdf <- "/gpfs01/home/shxyt1/data/sticklebacks/pheno/element_body_bubble_heatmap_sorted_by_height_rho.pdf"

# =========================================
# 2. User settings
# =========================================
missing_threshold <- 0.5

# Body trait columns in your table
body_traits <- c("body_length", "body_weight", "wet_weight_body_no_head")

# Labels shown in the plot
body_labels <- c(
  body_length = "Height",
  body_weight = "Weight",
  wet_weight_body_no_head = "Weight without head"
)

# Choose one body trait to define the x-axis order
# Options: "Height", "Weight", "Weight without head"
sort_trait_label <- "Height"

element_traits <- c(
  "Na","Mg","P","S","K","Ca","Al","V","Mn","Fe",
  "Co","Cu","Zn","As","Se","Rb","Sr","Mo","Ag","Cd","Cs","Ba","Tl","Pb","U"
)


# =========================================
# 3. Read data
# =========================================
dt <- fread(infile, sep = "\t", header = TRUE, data.table = FALSE)

element_traits <- intersect(element_traits, colnames(dt))
body_traits    <- intersect(body_traits, colnames(dt))

if (length(body_traits) == 0) {
  stop("None of the body trait columns were found. Please check body_traits names.")
}

for (v in c(element_traits, body_traits)) {
  dt[[v]] <- suppressWarnings(as.numeric(dt[[v]]))
}

# Keep only labels for body traits that exist
body_labels <- body_labels[names(body_labels) %in% body_traits]

if (!(sort_trait_label %in% body_labels)) {
  stop("sort_trait_label is not available in the current body trait labels.")
}

# =========================================
# 4. Filter elements by missingness
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

cat("Traits retained after missing-value filtering:\n")
print(element_traits_keep)

# =========================================
# 5. Build long element table and remove negatives
# =========================================
long_elem <- melt(
  as.data.table(dt[, element_traits_keep, drop = FALSE]),
  measure.vars = element_traits_keep,
  variable.name = "element",
  value.name = "value"
)

# Remove NA, non-finite values, and negative values
long_elem <- long_elem[!is.na(value) & is.finite(value) & value >= 0]

# Remove elements with too few values after filtering
trait_qc_after <- long_elem[, .(
  n_after_filter = .N,
  unique_after_filter = uniqueN(value)
), by = element]

elements_final <- trait_qc_after[
  n_after_filter > 3 & unique_after_filter > 1,
  element
]

if (length(elements_final) == 0) {
  stop("No element traits remained after filtering.")
}

cat("Traits retained after removing negatives:\n")
print(sort(elements_final))

# =========================================
# 6. Correlation function
# =========================================
get_cor_p <- function(x, y, method = "spearman") {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  
  if (length(x) < 4 || length(unique(x)) < 2 || length(unique(y)) < 2) {
    return(data.table(r = NA_real_, p = NA_real_, n = length(x)))
  }
  
  test <- suppressWarnings(cor.test(x, y, method = method, exact = FALSE))
  
  data.table(
    r = unname(test$estimate),
    p = test$p.value,
    n = length(x)
  )
}

# =========================================
# 7. Compute element-body correlations
# =========================================
res_list <- list()

for (bt in body_traits) {
  for (elem in elements_final) {
    tmp <- get_cor_p(dt[[bt]], dt[[elem]], method = "spearman")
    tmp[, body_trait := bt]
    tmp[, element := elem]
    res_list[[paste(bt, elem, sep = "__")]] <- tmp
  }
}

res_dt <- rbindlist(res_list)

# Add plotting labels
res_dt[, body_label := body_labels[body_trait]]

# =========================================
# 8. Order elements by one chosen body trait
#    Left to right:
#    strong positive -> weak positive -> weak negative -> strong negative
#    This is simply sorting rho from high to low
# =========================================
order_dt <- res_dt[body_label == sort_trait_label, .(element, r)]

# Put NA correlations at the end
order_dt <- order_dt[order(-r, element, na.last = TRUE)]

element_order <- order_dt$element

cat("Element order based on ", sort_trait_label, " rho:\n", sep = "")
print(order_dt)

# Apply factor levels
res_dt[, element := factor(element, levels = element_order)]

# Fix y-axis order
res_dt[, body_label := factor(
  body_label,
  levels = c("Weight without head", "Weight", "Height")
)]

# Keep only significant points for plotting
plot_dt <- res_dt[!is.na(r) & p < 0.05]

# =========================================
# 9. Plot
# =========================================
p <- ggplot(res_dt, aes(x = element, y = body_label)) +
  geom_tile(fill = "white", color = "#d9d9d9", linewidth = 0.5) +
  geom_point(
    data = plot_dt,
    aes(size = abs(r), color = r),
    shape = 16,
    alpha = 0.95
  ) +
  scale_size_continuous(
    range = c(0.4, 8),
    limits = c(0, 1),
    guide = "none"
  ) +
  scale_color_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Spearman rho"
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,    
      vjust = 1,    
      size = 11,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 13,
      color = "black",
      # face = "bold"
    ),
    legend.position = "right"
  )

print(p)

# =========================================
# 10. Save
# =========================================
plot_width <- max(10, 0.35 * length(element_order) + 5)


ggsave(out_png, p, width = plot_width, height = 2, dpi = 500, bg = "white")
ggsave(out_pdf, p, width = plot_width, height = 4.2, bg = "white")

cat("Saved to:\n")
cat(out_png, "\n")
cat(out_pdf, "\n")

