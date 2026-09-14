## Step 0 — Paths & parameters ================================================
library(data.table)

# Use cairo for bitmap devices (base R), avoid X11 font issues
options(bitmapType = "cairo")

# Input
infile <- "/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_nine_main/Six_log_main_TopSNPs_with_genes.csv"

# Output directory (requested)
outdir <- "/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_nine_main"

# Output files
gene_png  <- file.path(outdir, "SIGNIFICANT_gene_by_element_logp.png")
locus_png <- file.path(outdir, "SIGNIFICANT_locus250kbmerge_gene_by_element_logp.png")

# Locus rule: lead SNP distance <= 500 kb => same locus
# Implemented as merging +/- 250 kb windows
window_bp <- 250000L


## Step 1 — Load data & keep significant hits ================================
dt <- fread(infile)

# Keep significant hits only and require gene annotation
sig <- dt[Category == "significant" & !is.na(gene_symbol)]

# Ensure consistent types
sig[, `:=`(
  Element     = as.character(Element),
  gene_symbol = as.character(gene_symbol),
  CHR         = as.character(CHR),
  BP          = as.integer(BP),
  P           = as.numeric(P)
)]

# 1) List all elements that have at least one significant locus/gene
elements_with_sig <- sort(unique(sig$Element))
cat("Elements with significant loci:\n")
print(elements_with_sig)
cat("Number of elements with significant loci:", length(elements_with_sig), "\n")


## Step 2 — Heatmap saver (X=Element, Y=Gene; blue-white-red; legend) =========
save_heatmap_bwred <- function(mat, y_labels, x_labels, title, outfile,
                               base_cex = 1.8, axis_cex = 1.6, legend_cex = 1.2,
                               width = 2800, height = 1800) {
  stopifnot(is.matrix(mat))
  
  # Blue-White-Red palette
  pal <- colorRampPalette(c("blue", "white", "red"))(200)
  
  # Scale by observed max value (-log10 P), ignore NA
  zmax <- max(mat, na.rm = TRUE)
  if (!is.finite(zmax) || zmax <= 0) zmax <- 1
  breaks <- seq(0, zmax, length.out = length(pal) + 1)
  
  # Save with base R png using cairo backend
  png(outfile, width = width, height = height, type = "cairo")
  
  op <- par(no.readonly = TRUE)
  on.exit({ par(op); dev.off() }, add = TRUE)
  
  par(cex = base_cex)
  par(mar = c(7, 10, 4, 6))  # bottom, left, top, right
  
  # Make the first row appear at the TOP
  mat_plot <- mat[nrow(mat):1, , drop = FALSE]
  
  image(
    mat_plot,
    col = pal,
    breaks = breaks,
    axes = FALSE,
    xlab = "", ylab = "",
    main = title
  )
  
  # X axis: elements (horizontal labels)
  axis(1,
       at = seq(0, 1, length.out = length(x_labels)),
       labels = x_labels,
       las = 1, cex.axis = axis_cex)
  
  # Y axis: labels reversed to match mat_plot
  axis(2,
       at = seq(0, 1, length.out = length(y_labels)),
       labels = rev(y_labels),
       las = 2, cex.axis = axis_cex)
  
  box()
  
  # Legend color bar at top-right
  usr <- par("usr")
  x0 <- usr[2] + 0.02 * (usr[2] - usr[1])
  x1 <- usr[2] + 0.06 * (usr[2] - usr[1])
  y0 <- usr[4] - 0.35 * (usr[4] - usr[3])
  y1 <- usr[4] - 0.02 * (usr[4] - usr[3])
  
  yy <- seq(y0, y1, length.out = length(pal) + 1)
  for (i in 1:length(pal)) rect(x0, yy[i], x1, yy[i + 1], col = pal[i], border = NA, xpd = TRUE)
  rect(x0, y0, x1, y1, border = "black", xpd = TRUE)
  
  text(x1, y1, labels = "-log10(minP)", pos = 4, xpd = TRUE, cex = legend_cex)
  text(x1, y0, labels = sprintf("0 .. %.1f", zmax), pos = 4, xpd = TRUE, cex = legend_cex)
}




## Step 3 — Final polish (two-panel layout, 2-col legend, mild bold) ============

gene_png <- file.path(outdir, "genes_overlap_elements.png")

# --- (A) Build matrix: -log10(minP) per (gene, element) 
ge <- sig[, .(minP = min(P, na.rm = TRUE)), by = .(gene_symbol, Element)]
ge[, logp := -log10(minP)]
gene_dt <- dcast(ge, gene_symbol ~ Element, value.var = "logp", fill = NA_real_)

genes <- gene_dt$gene_symbol
elements <- setdiff(colnames(gene_dt), "gene_symbol")
mat <- as.matrix(gene_dt[, ..elements])

# --- (B) Order genes by overlap, cluster within overlap group 
B <- !is.na(mat)
overlap_n <- rowSums(B)

groups <- split(seq_len(nrow(mat)), overlap_n)
grp_keys <- sort(as.integer(names(groups)), decreasing = TRUE)
ord_all <- integer(0)
for (k in grp_keys) {
  idx <- groups[[as.character(k)]]
  if (length(idx) <= 2) ord_all <- c(ord_all, idx) else {
    D <- dist(B[idx, , drop = FALSE], method = "binary")
    hc <- hclust(D, method = "average")
    ord_all <- c(ord_all, idx[hc$order])
  }
}
mat <- mat[ord_all, , drop = FALSE]
genes <- genes[ord_all]

# Order elements by number of significant genes (optional)
elem_hits <- colSums(!is.na(mat))
elem_ord <- order(elem_hits, decreasing = TRUE)
mat <- mat[, elem_ord, drop = FALSE]
elements <- elements[elem_ord]

# --- (C) Bins start near genome-wide threshold: p < 5e-8 => -log10(p) > 7.30 ---
# Use neat bins starting from 7.5
cuts <- c(7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, Inf)  # 8 bins

bin_id <- matrix(0L, nrow(mat), ncol(mat))  # 0=absent
for (i in seq_len(nrow(mat))) for (j in seq_len(ncol(mat))) {
  if (!is.na(mat[i, j])) {
    v <- mat[i, j]
    # clamp anything below first cut into bin 1 (still lightest)
    bin_id[i, j] <- if (v <= cuts[1]) 1L else max(1L, min(8L, findInterval(v, cuts, rightmost.closed = TRUE)))
  }
}

# --- (D) Colors: 8 orange->red shades; absent=white 
anchors <- c("#fee0b6", "#fdae6b", "#e6550d", "#a50f15")
reds8 <- colorRampPalette(anchors)(8)

cols <- c("white", reds8)

legend_title <- "-log10(P)"
legend_labels <- c("(7.5,8.0]", "(8.0,8.5]", "(8.5,9.0]", "(9.0,9.5]",
                   "(9.5,10.0]", "(10.0,10.5]", "(10.5,11.0]", "(11.0,+)")

# --- (E) Device: stable size; do NOT overblow fonts 
png(gene_png, width = 4200, height = max(3200, 55 * length(genes)),
    res = 300, pointsize = 12, type = "cairo")

op <- par(no.readonly = TRUE)
on.exit({ par(op); dev.off() }, add = TRUE)

names(grDevices::quartzFonts())  # macOS 才有

par(family = "Arial")  # will fall back if Arial not available

# --- (F) Layout: heatmap + legend panel (legend in 2 columns) 
layout(matrix(c(1,2), nrow = 1), widths = c(4.4, 1.2))
par(oma = c(0, 0, 2.2, 0))  # outer top margin for title

nr <- nrow(bin_id); nc <- ncol(bin_id)

# Font sizes: keep ticks as baseline; titles only slightly bigger
tick_cex    <- min(1.4, 36 / max(20, nr)) * 2.4    # baseline (tick labels)
axislab_cex <- tick_cex * 1.15                     # slightly bigger than ticks
title_cex   <- tick_cex * 1.75                     # slightly bigger than ticks
leg_cex     <- tick_cex                             # legend same as ticks

# --- Panel 1: heatmap 
par(mar = c(5.5, 18, 1.0, 0.8), bty = "n", xaxs = "i", yaxs = "i", xpd = NA)

Z <- bin_id[nr:1, , drop = FALSE]
x <- seq_len(nc); y <- seq_len(nr)

image(x = x, y = y, z = t(Z),
      col = cols, breaks = seq(-0.5, 8.5, by = 1),
      axes = FALSE, xlab = "", ylab = "", main = "")

# Title (mild bold)
mtext("Significant genes across elements", outer = TRUE, side = 3, line = 0.6,
      cex = title_cex, font = 2)

# Axis titles: mild bold, but not huge
mtext("Element", side = 1, line = 3.2, cex = axislab_cex, font = 2)
# Push "Gene" left to avoid overlapping gene labels
mtext("Gene", side = 2, line = 12.5, cex = axislab_cex, font = 2)


# Tick labels: normal weight (not bold)
axis(1, at = x, labels = elements, las = 1, cex.axis = tick_cex, lwd = 0, lwd.ticks = 0, font = 1)
axis(2, at = y, labels = rev(genes), las = 2, cex.axis = tick_cex, lwd = 0, lwd.ticks = 0, font = 1)

# --- Panel 2: legend (2 columns) -
par(mar = c(5.5, 0.5, 1.0, 1.0), bty = "n", xaxs = "i", yaxs = "i")
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))

# Legend title: mild bold, similar size
text(0.05, 0.95, legend_title, adj = 0, cex = leg_cex, font = 2)

# Two-column legend layout
# We have 9 items (absent + 8 bins). Put 5 in col1, 4 in col2.
items <- 2:(length(cols))
col1 <- items[1:4]
col2 <- items[5:8]


sw_w <- 0.045
row_step <- 0.035

y_top <- 0.86
x1 <- 0.05
x2 <- 0.55

draw_leg <- function(idxs, x0) {
  for (ii in seq_along(idxs)) {
    k <- idxs[ii]
    yy <- y_top - (ii - 1) * row_step
    symbols(x0 + sw_w/2, yy,
            squares = sw_w,
            inches = FALSE,
            add = TRUE,
            bg = cols[k],
            fg = NA)
    text(x0 + sw_w + 0.05, yy, legend_labels[k - 1], adj = 0, cex = leg_cex, font = 1)
  }
}
draw_leg(col1, x1)
draw_leg(col2, x2)

dev.off()
cat("Saved:", gene_png, "\n")


## Step 4 — Locus-level overlap (±250kb; lead SNP distance <= 500kb) ============
library(data.table)

setDT(sig)

# Make sure types are correct
sig[, BP := as.integer(BP)]
sig[, CHR := as.character(CHR)]

# 1) Diagnose missing values
cat("NA in BP:", sum(is.na(sig$BP)), "\n")
cat("NA in CHR:", sum(is.na(sig$CHR)), "\n")
cat("NA in Element:", sum(is.na(sig$Element)), "\n")

# 2) Drop rows with missing key fields (required for locus building)
sig <- sig[!is.na(Element) & !is.na(CHR) & !is.na(BP)]

# 3) Sort (this is REQUIRED)
setorder(sig, Element, CHR, BP)






#Step 4.1Within each element, sort and group SNPS by position (≤500kb)
library(data.table)

setDT(sig)

# Ensure numeric
sig[, BP := as.integer(BP)]
sig[, CHR := as.character(CHR)]

# Sort by element, chromosome, position
setorder(sig, Element, CHR, BP)

# Assign locus within each element + chromosome
sig[, locus_tmp := {
  d <- c(Inf, diff(BP))          # first SNP always starts a new locus
  cumsum(d > 500000) + 1L        # locus index starts at 1
}, by = .(Element, CHR)]


#Step 4.2 - Summarize into a locus level table
locus_table <- sig[, .(
  start        = min(BP),
  end          = max(BP),
  lead_snp     = SNP[which.min(P)],
  bestP        = min(P),
  best_gene    = gene_symbol[which.min(P)],
  genes_in_locus = paste(sort(unique(gene_symbol)), collapse = ",")
), by = .(Element, CHR, locus_tmp)]
t

#Step 4.3 —— Number the locus (within each element)

setorder(locus_table, Element, CHR, start)

locus_table[, locus_id := seq_len(.N), by = Element]

# Reorder columns
setcolorder(
  locus_table,
  c("Element", "locus_id", "CHR", "start", "end",
    "lead_snp", "bestP", "best_gene", "genes_in_locus")
)

#Step 4.4: Save to your specified directory
outdir <- "/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_nine_main"
out_locus <- file.path(outdir, "locus_table_500kb.tsv")

fwrite(locus_table, out_locus, sep = "\t")
cat("Saved locus table to:", out_locus, "\n")



#Step 5locus-level overlap heatmap=====
library(data.table)

outdir <- "/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_nine_main"

locus <- fread(file.path(outdir, "locus_table_500kb.tsv"))
locus[, logp := -log10(bestP)]

# Expand genes_in_locus into long format: one row per (locus, gene)
lg <- locus[, .(gene_symbol = unique(unlist(strsplit(genes_in_locus, "\\s*,\\s*")))),
            by = .(Element, CHR, start, end, locus_id, bestP, logp, lead_snp, best_gene)]

# Drop empty gene strings (safety)
lg <- lg[!is.na(gene_symbol) & gene_symbol != ""]

#Step 5.1 - Construct the locus × element matrix

ge <- lg[, .(
  minP = min(bestP, na.rm = TRUE),
  logp = max(logp,  na.rm = TRUE)   # equivalent to -log10(minP)
), by = .(gene_symbol, Element)]

gene_dt <- dcast(ge, gene_symbol ~ Element, value.var = "logp", fill = NA_real_)

genes <- gene_dt$gene_symbol
elements <- setdiff(colnames(gene_dt), "gene_symbol")
mat <- as.matrix(gene_dt[, ..elements])
rownames(mat) <- genes

#Step 5.2 - Sort by "whether hit by multiple elements" 
# Order loci by maximum signal strength across elements
B <- !is.na(mat)
overlap_n <- rowSums(B)

groups <- split(seq_len(nrow(mat)), overlap_n)
grp_keys <- sort(as.integer(names(groups)), decreasing = TRUE)

ord_all <- integer(0)
for (k in grp_keys) {
  idx <- groups[[as.character(k)]]
  if (length(idx) <= 2) {
    ord_all <- c(ord_all, idx)
  } else {
    D <- dist(B[idx, , drop = FALSE], method = "binary")
    hc <- hclust(D, method = "average")
    ord_all <- c(ord_all, idx[hc$order])
  }
}

mat <- mat[ord_all, , drop = FALSE]
genes <- genes[ord_all]

# Order elements by number of hit genes (optional)
elem_hits <- colSums(!is.na(mat))
elem_ord <- order(elem_hits, decreasing = TRUE)
mat <- mat[, elem_ord, drop = FALSE]
elements <- elements[elem_ord]

#Step5.3 picture
gene_png <- file.path(outdir, "genes_overlap_elements.png")

# ---- bins: start from genome-wide threshold -log10(5e-8)=7.301... -------------
cuts <- c(7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11.0, Inf)  # 8 bins

bin_id <- matrix(0L, nrow(mat), ncol(mat))  # 0=absent
for (i in seq_len(nrow(mat))) {
  for (j in seq_len(ncol(mat))) {
    if (!is.na(mat[i, j])) {
      v <- mat[i, j]
      # clamp values below first cut into lightest bin
      if (v <= cuts[1]) {
        bin_id[i, j] <- 1L
      } else {
        bin_id[i, j] <- min(8L, findInterval(v, cuts, rightmost.closed = TRUE))
      }
    }
  }
}

# ---- colors: 8 orange->red shades --------------------------------------------
anchors <- c("#fee0b6", "#fdae6b", "#e6550d", "#a50f15")
reds8 <- colorRampPalette(anchors)(8)
cols <- c("white", reds8)

legend_title <- "-log10(P)"
legend_labels <- c("(7.5,8.0]", "(8.0,8.5]", "(8.5,9.0]", "(9.0,9.5]",
                   "(9.5,10.0]", "(10.0,10.5]", "(10.5,11.0]", "(11.0,+)")

# ---- device ------------------------------------------------------------------
png(gene_png,
    width = 4200,
    height = max(3200, 55 * length(genes)),
    res = 300,
    type = "cairo")

op <- par(no.readonly = TRUE)
on.exit({ par(op); dev.off() }, add = TRUE)

# Font family: use Arial if available; otherwise fallback to sans
par(family = "sans")

# ---- layout: heatmap + legend ------------------------------------------------
layout(matrix(c(1,2), nrow = 1), widths = c(4.4, 1.2))
par(oma = c(0, 0, 2.0, 0))

nr <- nrow(bin_id); nc <- ncol(bin_id)

# ---- font sizes: YOU adjust here ---------------------------------------------
tick_cex    <- 1.15   # tick labels (genes + elements)
axislab_cex <- 1.30   # axis titles
title_cex   <- 1.45   # main title
leg_cex     <- 1.15   # legend text

# ---- panel 1: heatmap --------------------------------------------------------
par(mar = c(5.5, 14.5, 1.0, 0.8), bty = "n", xaxs = "i", yaxs = "i", xpd = NA)

Z <- bin_id[nr:1, , drop = FALSE]
x <- seq_len(nc); y <- seq_len(nr)

image(x = x, y = y, z = t(Z),
      col = cols, breaks = seq(-0.5, 8.5, by = 1),
      axes = FALSE, xlab = "", ylab = "", main = "")

mtext("Significant genes across elements", outer = TRUE, side = 3, line = 0.5,
      cex = title_cex, font = 2)

mtext("Element", side = 1, line = 3.2, cex = axislab_cex, font = 2)
mtext("Gene",    side = 2, line = 9.5, cex = axislab_cex, font = 2)

axis(1, at = x, labels = elements, las = 1, cex.axis = tick_cex, lwd = 0, lwd.ticks = 0)
axis(2, at = y, labels = rev(genes), las = 2, cex.axis = tick_cex, lwd = 0, lwd.ticks = 0)

# ---- panel 2: legend (2 columns, 4 + 4; no "absent") -------------------------
par(mar = c(5.5, 0.5, 1.0, 1.0), bty = "n", xaxs = "i", yaxs = "i")
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))

text(0.05, 0.95, legend_title, adj = 0, cex = leg_cex, font = 2)

# exactly 8 bins -> 4 per column
items <- 2:(length(cols))  # bins only, exclude "white"
col1 <- items[1:4]
col2 <- items[5:8]

# make true squares in user coords: use symbols(..., squares=..., inches=TRUE)
# inches=TRUE makes squares in inches (device space), always square.
sw_in <- 0.12
row_step <- 0.07
y_top <- 0.83
x1 <- 0.05
x2 <- 0.55

draw_leg <- function(idxs, x0, lab_idx_start) {
  for (ii in seq_along(idxs)) {
    k <- idxs[ii]
    yy <- y_top - (ii - 1) * row_step
    symbols(x0, yy, squares = sw_in, inches = TRUE, add = TRUE, bg = cols[k], fg = NA)
    text(x0 + 0.09, yy, legend_labels[lab_idx_start + ii - 1], adj = 0, cex = leg_cex, font = 1)
  }
}
draw_leg(col1, x1, 1)
draw_leg(col2, x2, 5)

dev.off()
cat("Saved:", gene_png, "\n")




