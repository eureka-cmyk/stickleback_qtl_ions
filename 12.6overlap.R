library(data.table)
options(bitmapType = "cairo")

# 1. Load GWAS top SNPs with gene annotations
dt <- fread("/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_no_outlier/All_log_TopSNPs_with_genes.csv")

# 2. Keep only rows with valid gene annotations
dt <- dt[!is.na(gene_symbol)]

# 3. Focus on genome-wide significant signals only
dt_sig <- dt[Category == "significant"]

# 4. Create a list of significant genes for each element
elem_gene_list <- split(unique(dt_sig$gene_symbol), dt_sig$Element)

# 5. Count how many elements each gene is significant in
gene_elem_count <- dt_sig[, .(
  n_element = uniqueN(Element),
  elements  = paste(sort(unique(Element)), collapse = ",")
), by = gene_symbol]

# 6. Extract overlapping genes (significant in >= 2 elements)
overlap_genes <- gene_elem_count[n_element >= 2][order(-n_element)]

# Inspect results
print(overlap_genes)


# Keep only gene symbols that overlap across >=2 elements
overlap_gene_set <- overlap_genes$gene_symbol

# Build gene x element presence/absence matrix (overlap genes only)
ge <- unique(dt_sig[gene_symbol %in% overlap_gene_set, .(gene_symbol, Element)])
ge[, present := 1L]
gene_mat_dt <- dcast(ge, gene_symbol ~ Element, value.var = "present", fill = 0L)

gene_names <- gene_mat_dt$gene_symbol
gene_mat <- as.matrix(gene_mat_dt[, -"gene_symbol"])

# Order genes by number of elements (descending)
ord_g <- order(rowSums(gene_mat), decreasing = TRUE)
gene_mat <- gene_mat[ord_g, , drop = FALSE]
gene_names <- gene_names[ord_g]

# Plot a simple presence/absence heatmap (base R)
gene_png <- "overlap_significant_gene_level.png"
png(gene_png, width = 1600, height = max(900, 25 * nrow(gene_mat)), type = "cairo")
par(mar = c(10, 20, 4, 2))
image(
  t(gene_mat[nrow(gene_mat):1, , drop = FALSE]),
  axes = FALSE,
  xlab = "", ylab = "",
  main = "Overlap of Significant Hits (Gene-level)"
)
par(family = "sans")

axis(1, at = seq(0, 1, length.out = nrow(gene_mat)), labels = rev(gene_names), las = 2, cex.axis = 0.6)
axis(2, at = seq(0, 1, length.out = ncol(gene_mat)), labels = colnames(gene_mat), las = 2, cex.axis = 0.9)
box()
dev.off()
cat("Saved gene-level plot:", gene_png, "\n")


# ============================================================
# B) Locus-level overlap using ±50kb windows (significant only)
#    1) Expand each SNP into a window
#    2) Merge overlapping windows per chromosome into loci
#    3) Compute locus overlap across elements and plot heatmap
# ============================================================

window_bp <- 50000L

# Create +/-50kb windows around each significant SNP
w <- dt_sig[, .(
  CHR = as.character(CHR),
  start = pmax(1L, as.integer(BP) - window_bp),
  end   = as.integer(BP) + window_bp,
  Element,
  SNP
)]
setorder(w, CHR, start, end)

# Helper: merge overlapping windows within each chromosome
merge_chr_windows <- function(sub) {
  # sub must be sorted by start
  loci <- vector("list", 0L)
  
  cur_start <- sub$start[1]
  cur_end   <- sub$end[1]
  cur_elems <- unique(sub$Element[1])
  cur_snps  <- unique(sub$SNP[1])
  
  if (nrow(sub) > 1) {
    for (i in 2:nrow(sub)) {
      s <- sub$start[i]
      e <- sub$end[i]
      el <- sub$Element[i]
      sn <- sub$SNP[i]
      
      if (s <= cur_end) {
        # Overlap: extend current locus
        cur_end <- max(cur_end, e)
        cur_elems <- unique(c(cur_elems, el))
        cur_snps  <- unique(c(cur_snps, sn))
      } else {
        # No overlap: finalize current locus
        loci[[length(loci) + 1L]] <- list(
          start = cur_start,
          end   = cur_end,
          elements = paste(sort(cur_elems), collapse = ","),
          n_element = length(cur_elems),
          snps = paste(cur_snps, collapse = ";")
        )
        # Start a new locus
        cur_start <- s
        cur_end   <- e
        cur_elems <- unique(el)
        cur_snps  <- unique(sn)
      }
    }
  }
  
  # Add the last locus
  loci[[length(loci) + 1L]] <- list(
    start = cur_start,
    end   = cur_end,
    elements = paste(sort(cur_elems), collapse = ","),
    n_element = length(cur_elems),
    snps = paste(cur_snps, collapse = ";")
  )
  
  rbindlist(loci)
}

# Merge windows to loci by chromosome
locus_dt <- w[, merge_chr_windows(.SD), by = CHR]
locus_dt[, locus_id := paste0("chr", CHR, ":", start, "-", end)]

# Keep only overlapping loci (shared by >=2 elements)
locus_overlap <- locus_dt[n_element >= 2]
if (nrow(locus_overlap) == 0) stop("No overlapping loci found at ±50kb (shared by >=2 elements).")

# Build locus x element presence/absence matrix
exp <- locus_overlap[, .(Element = unlist(strsplit(elements, ",", fixed = TRUE))), by = locus_id]
exp[, present := 1L]
locus_mat_dt <- dcast(exp, locus_id ~ Element, value.var = "present", fill = 0L)

locus_names <- locus_mat_dt$locus_id
locus_mat <- as.matrix(locus_mat_dt[, -"locus_id"])

# Order loci by number of elements (descending)
ord_l <- order(rowSums(locus_mat), decreasing = TRUE)
locus_mat <- locus_mat[ord_l, , drop = FALSE]
locus_names <- locus_names[ord_l]

# Plot locus-level presence/absence heatmap
locus_png <- "overlap_significant_locus_level_50kb.png"
png(locus_png, width = 1800, height = max(900, 22 * nrow(locus_mat)))
par(mar = c(10, 22, 4, 2))
image(
  t(locus_mat[nrow(locus_mat):1, , drop = FALSE]),
  axes = FALSE,
  xlab = "", ylab = "",
  main = "Overlap of Significant Hits (Locus-level, ±50kb)"
)
axis(1, at = seq(0, 1, length.out = nrow(locus_mat)), labels = rev(locus_names), las = 2, cex.axis = 0.55)
axis(2, at = seq(0, 1, length.out = ncol(locus_mat)), labels = colnames(locus_mat), las = 2, cex.axis = 0.9)
box()
dev.off()
cat("Saved locus-level plot:", locus_png, "\n")

# Optional: print locus overlap summary table
locus_summary <- data.table(locus_id = locus_names, n_element = rowSums(locus_mat))
print(locus_summary[order(-n_element)])
