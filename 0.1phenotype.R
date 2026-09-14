library(data.table)
library(ggplot2)
library(patchwork)

options(bitmapType = "cairo")

# ================================
# 1. Input and output
# ================================
infile <- "/gpfs01/home/shxyt1/data/sticklebacks/pheno/all_pheno_clean.txt"
outpdf <- "/gpfs01/home/shxyt1/data/sticklebacks/pheno/all_pheno_raw_distributions_4x4.pdf"

# ================================
# 2. Read data
# ================================
dt <- fread(infile, sep = "\t", header = TRUE, data.table = FALSE)

# ================================
# 3. Columns to exclude
# ================================
exclude_cols <- c(
  "FID", "IID",
  "date_fish_caught", "date_fish_processed",
  "sex", "schisto", "notes", "icpms_run",
  "Plates", "1st spine length", "2nd spine length", "Pelvis spine length", "age"
)

# Keep only columns that actually exist
exclude_cols <- intersect(exclude_cols, colnames(dt))

# Candidate columns
candidate_cols <- setdiff(colnames(dt), exclude_cols)

# ================================
# 4. Force candidate columns to numeric if possible
# ================================
for (col in candidate_cols) {
  dt[[col]] <- suppressWarnings(as.numeric(dt[[col]]))
}

# Keep only numeric columns
numeric_cols <- candidate_cols[sapply(dt[candidate_cols], is.numeric)]

# Remove columns that are all NA or have <=1 unique non-NA value
numeric_cols <- numeric_cols[
  sapply(numeric_cols, function(col) {
    x <- dt[[col]]
    x <- x[!is.na(x)]
    length(x) > 1 && length(unique(x)) > 1
  })
]

cat("Number of numeric variables to plot:", length(numeric_cols), "\n")
cat("Variables to plot:\n")
print(numeric_cols)

# ================================
# 5. Function to make one distribution plot
# ================================
make_dist_plot <- function(data, varname) {
  df <- data.frame(value = data[[varname]])
  df <- df[!is.na(df$value), , drop = FALSE]
  
  p <- ggplot(df, aes(x = value)) +
    geom_histogram(aes(y = after_stat(density)),
                   bins = 30,
                   fill = "#6BAED6",
                   color = "black",
                   linewidth = 0.2) +
    geom_density(color = "black", linewidth = 0.7, na.rm = TRUE) +
    labs(
      title = varname,
      x = NULL,
      y = "Density"
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10, face = "plain"),
      axis.title.y = element_text(size = 9),
      axis.text = element_text(size = 8, color = "black"),
      axis.line.x = element_line(color = "black", linewidth = 0.6),
      axis.line.y = element_line(color = "black", linewidth = 0.6),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      plot.margin = margin(4, 4, 4, 4)
    )
  
  return(p)
}

# ================================
# 6. Build all plots
# ================================
plot_list <- lapply(numeric_cols, function(col) make_dist_plot(dt, col))

# ================================
# 7. Save to PDF in 4 x 4 layout
# ================================
plots_per_page <- 16
n_plots <- length(plot_list)
n_pages <- ceiling(n_plots / plots_per_page)

pdf(outpdf, width = 14, height = 14, onefile = TRUE)

for (i in seq_len(n_pages)) {
  start_idx <- (i - 1) * plots_per_page + 1
  end_idx <- min(i * plots_per_page, n_plots)
  
  subplots <- plot_list[start_idx:end_idx]
  
  # If last page has fewer than 16 plots, fill with blank plots
  if (length(subplots) < plots_per_page) {
    n_missing <- plots_per_page - length(subplots)
    blank_plots <- replicate(
      n_missing,
      ggplot() + theme_void(),
      simplify = FALSE
    )
    subplots <- c(subplots, blank_plots)
  }
  
  combined_plot <- wrap_plots(subplots, ncol = 4, nrow = 4)
  print(combined_plot)
}

dev.off()

cat("PDF saved to:\n", outpdf, "\n")

