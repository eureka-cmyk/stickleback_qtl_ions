# ============================================================
# Generate log(LOD-replaced) phenotype files + distributions + PCA
# Author: (your name)
# Notes:
#   - Replace negative element values with element-specific LOD
#   - Then apply natural log() transform
#   - Output per-element phenotype files: FID IID PHENO
#   - Make distribution plots (9 plots per page)
#   - Run PCA on all log_LOD elements and export PC1-5
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(ggplot2)
  library(patchwork)
  
})

# -----------------------------
# Paths (edit only if needed)
# -----------------------------
lod_xlsx <- "/gpfs01/home/shxyt1/data/sticklebacks/pheno/Limits of detection.xlsx"
pheno_in <- "/gpfs01/home/shxyt1/data/sticklebacks/pheno/all_pheno_clean.txt"
out_dir  <- "/gpfs01/home/shxyt1/data/sticklebacks/pheno/traits/log_LOD"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Read main phenotype table
# -----------------------------
dt <- fread(pheno_in, sep = "\t", header = TRUE, data.table = TRUE)

# Ensure FID/IID exist
stopifnot(all(c("FID","IID") %in% names(dt)))

# -----------------------------
# Read LOD table from Excel
# Strategy:
#   - Read entire sheet as text-ish
#   - Find the row containing "Limit of detection"
#   - Extract numeric LOD values across element columns
# -----------------------------
lod_raw <- read_excel(lod_xlsx, sheet = 1, col_names = TRUE)

lod_dt <- as.data.table(lod_raw)

# Find a row that contains the keyword anywhere
keyword <- "Limit of detection"
row_has_keyword <- apply(lod_dt, 1, function(x) any(grepl(keyword, as.character(x), ignore.case = TRUE)))
if (!any(row_has_keyword)) {
  stop("Cannot find a row containing 'Limit of detection' in the Excel file. Please check the sheet.")
}
lod_row_idx <- which(row_has_keyword)[1]

# Identify candidate element columns (intersection with dt columns)
# Keep only columns that look like element symbols present in dt
element_cols <- intersect(names(dt), names(lod_dt))

if (length(element_cols) < 5) {
  stop("Too few overlapping element columns between all_pheno_clean.txt and the LOD Excel sheet. Check column names.")
}

# Extract the LOD row for those element columns
lod_vals <- lod_dt[lod_row_idx, ..element_cols]

# Convert to numeric (handle character cells)
lod_numeric <- as.list(lod_vals)
lod_numeric <- lapply(lod_numeric, function(x) suppressWarnings(as.numeric(x)))
lod_numeric <- unlist(lod_numeric, use.names = TRUE)

# Remove elements with missing LOD values
missing_lod <- names(lod_numeric)[is.na(lod_numeric)]
if (length(missing_lod) > 0) {
  message("Warning: Missing LOD for these elements (will be skipped): ", paste(missing_lod, collapse = ", "))
  lod_numeric <- lod_numeric[!is.na(lod_numeric)]
}

elements <- names(lod_numeric)

# -----------------------------
# Create LOD-replaced dataset (before log)
# -----------------------------
dt_lod <- copy(dt)

for (el in elements) {
  x <- dt_lod[[el]]
  
  # Replace negative values with element-specific LOD (before log)
  x[!is.na(x) & x < 0] <- lod_numeric[[el]]
  
  dt_lod[[el]] <- x
}

message("Elements with LOD found: ", length(elements))
message("First few elements: ", paste(head(elements, 10), collapse = ", "))

# -----------------------------
# Replace negatives with LOD and log-transform
# Output per-element phenotype files
# -----------------------------
replaced_summary <- data.table(
  element = elements,
  n_total = NA_integer_,
  n_negative = NA_integer_,
  n_replaced = NA_integer_
)

for (el in elements) {
  x <- dt[[el]]
  
  # Count negatives (ignore NA)
  n_total    <- sum(!is.na(x))
  n_negative <- sum(!is.na(x) & x < 0)
  
  # Replace negatives with LOD
  x2 <- x
  x2[!is.na(x2) & x2 < 0] <- lod_numeric[[el]]
  
  n_replaced <- sum(!is.na(x) & x < 0)
  
  # Natural log transform
  # Note: after replacement, values should be >0 (LOD is positive)
  xlog <- rep(NA_real_, length(x2))
  ok <- !is.na(x2) & x2 > 0
  xlog[ok] <- log(x2[ok])
  
  # Write phenotype file: FID IID PHENO
  out_file <- file.path(out_dir, paste0("pheno_", el, "_log_LOD.txt"))
  out_dt <- dt[, .(FID, IID)]
  out_dt[, PHENO := xlog]
  fwrite(out_dt, out_file,
         sep = "\t",
         quote = FALSE,
         na = "NA",
         col.names = FALSE)
  
  
  # Save summary
  replaced_summary[element == el, `:=`(n_total = n_total,
                                       n_negative = n_negative,
                                       n_replaced = n_replaced)]
}

# Save replacement summary
fwrite(replaced_summary, file.path(out_dir, "LOD_replacement_summary.txt"),
       sep = "\t", quote = FALSE)

message("Finished writing per-element log_LOD phenotype files to: ", out_dir)

# -----------------------------
# Distribution plots (9 per page)
# Output: log_LOD_distributions_9in1.pdf
# -----------------------------
# We will try to use patchwork if installed; otherwise fallback to gridExtra
has_patchwork <- requireNamespace("patchwork", quietly = TRUE)
has_gridExtra <- requireNamespace("gridExtra", quietly = TRUE)

plot_list <- list()

for (el in elements) {
  f <- file.path(out_dir, paste0("pheno_", el, "_log_LOD.txt"))
  d <- fread(f, sep = "\t", header = TRUE)
  p <- ggplot(d, aes(x = PHENO)) +
    geom_histogram(aes(y = after_stat(density)), bins = 40, fill = "#6BAED6", color = "black", alpha = 0.7) +
    geom_density(linewidth = 0.8) +
    labs(title = el, x = NULL, y = "Density") +
    theme_bw(base_size = 12)+
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"))
  plot_list[[el]] <- p
}

pdf_out <- file.path(out_dir, "log_LOD_distributions_9in1.pdf")
pdf(pdf_out, width = 11, height = 8.5)

if (has_patchwork) {
  # 9 plots per page with patchwork
  library(patchwork)
  els <- names(plot_list)
  pages <- split(els, ceiling(seq_along(els) / 9))
  for (pg in pages) {
    
    # get plots for this page
    p_list <- plot_list[pg]
    
    # add empty plots to make 9 panels
    n_missing <- 9 - length(p_list)
    if (n_missing > 0) {
      empty_plots <- replicate(n_missing, patchwork::plot_spacer(), simplify = FALSE)
      p_list <- c(p_list, empty_plots)
    }
    
    # fixed 3x3 layout
    p <- wrap_plots(p_list, ncol = 3)
    
    print(p)
  }
} else if (has_gridExtra) {
  # 9 plots per page with gridExtra
  library(gridExtra)
  els <- names(plot_list)
  pages <- split(els, ceiling(seq_along(els) / 9))
  for (pg in pages) {
    grid.arrange(grobs = plot_list[pg], ncol = 3)
  }
} else {
  # Fallback: one plot per page
  message("Neither patchwork nor gridExtra found; saving one plot per page.")
  for (el in names(plot_list)) print(plot_list[[el]])
}

dev.off()
message("Saved distribution PDF: ", pdf_out)

# -----------------------------
# PCA using all log_LOD elements
# Steps:
#   - build a matrix: rows=individuals, cols=elements
#   - keep complete cases (no NA)
#   - scale and run prcomp
#   - export PC1-5 with FID IID
# -----------------------------
# Build log_LOD matrix
log_mat <- dt[, ..elements]
# Replace negatives and log again (to ensure PCA consistent with outputs)
for (el in elements) {
  x <- log_mat[[el]]
  x[!is.na(x) & x < 0] <- lod_numeric[[el]]
  xlog <- rep(NA_real_, length(x))
  ok <- !is.na(x) & x > 0
  xlog[ok] <- log(x[ok])
  log_mat[[el]] <- xlog
}

# Keep rows with complete data for PCA
keep <- complete.cases(log_mat)
message("PCA complete cases: ", sum(keep), " / ", nrow(log_mat))

pca_res <- prcomp(as.matrix(log_mat[keep, ]), center = TRUE, scale. = TRUE)

pcs <- as.data.table(pca_res$x[, 1:5, drop = FALSE])
setnames(pcs, paste0("PC", 1:5))

pca_out <- dt[keep, .(FID, IID)]
pca_out <- cbind(pca_out, pcs)

pca_file <- file.path(out_dir, "log_LOD_PCA_PC1-5.txt")
fwrite(pca_out, pca_file,
       sep = "\t",
       quote = FALSE,
       na = "NA",
       col.names = FALSE)

message("Saved PCA PC1-5: ", pca_file)



# ============================================================
# Distribution plots AFTER LOD replacement (before log)
# ============================================================

# -----------------------------
# Save full phenotype table after LOD replacement (before log)
# -----------------------------
lod_table_out <- "/gpfs01/home/shxyt1/data/sticklebacks/pheno/LOD_replaced_pheno.txt"
fwrite(dt_lod, lod_table_out, sep = "\t", quote = FALSE, na = "NA")
message("Saved LOD-replaced phenotype table (before log): ", lod_table_out)

# -----------------------------
# Distribution plots AFTER LOD replacement (before log)
# Output: LOD_replaced_distributions_9in1.pdf
# -----------------------------
plot_list_lod <- list()

for (el in elements) {
  
  d <- dt_lod[, .(PHENO = get(el))]
  
  p <- ggplot(d, aes(x = PHENO)) +
    geom_histogram(aes(y = after_stat(density)),
                   bins = 40,
                   fill = "#6BAED6",
                   color = "black",
                   alpha = 0.7) +
    geom_density(linewidth = 0.8) +
    labs(title = el, x = NULL, y = "Density") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  plot_list_lod[[el]] <- p
}

pdf_out_lod <- file.path(out_dir, "LOD_replaced_distributions_9in1.pdf")
pdf(pdf_out_lod, width = 11, height = 8.5)

els <- names(plot_list_lod)
pages <- split(els, ceiling(seq_along(els) / 9))

for (pg in pages) {
  
  p_list <- plot_list_lod[pg]
  
  # Fill empty slots to keep identical panel sizes on the last page
  n_missing <- 9 - length(p_list)
  if (n_missing > 0) {
    empty_plots <- replicate(n_missing, patchwork::plot_spacer(), simplify = FALSE)
    p_list <- c(p_list, empty_plots)
  }
  
  print(wrap_plots(p_list, ncol = 3))
}

dev.off()
message("Saved LOD-replaced distributions (before log): ", pdf_out_lod)
