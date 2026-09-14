library(data.table)
options(bitmapType = "cairo")

# Read your P < 1e-6 SNP table
dt <- fread("/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_nine_main/Six_log_SNPs_with_genes_1e-6.csv")

# Ensure correct types
dt[, CHR := as.integer(CHR)]
dt[, BP  := as.integer(BP)]
dt[, P   := as.numeric(P)]

# -----------------------------
# Parameters you may change
# -----------------------------
p_cutoff <- 1e-6
clump_kb <- 500  # set to 250 or 500
clump_bp <- clump_kb * 1000

# Keep only SNPs passing the threshold (in case the input file contains others)
dt <- dt[P <= p_cutoff]

# Sort within each Element × CHR by significance first
setorder(dt, Element, CHR, P, BP)

# -----------------------------
# Distance-based clumping:
# pick the most significant SNP as lead,
# assign all SNPs within ±clump_bp to this locus,
# repeat on remaining SNPs.
# -----------------------------
get_distance_loci <- function(d) {
  d <- copy(d)
  d[, assigned := FALSE]
  locus_idx <- 0L
  
  while (any(d$assigned == FALSE)) {
    locus_idx <- locus_idx + 1L
    
    # Pick the lead SNP among unassigned SNPs (smallest P-value)
    lead_row <- d[assigned == FALSE][which.min(P)]
    lead_bp  <- lead_row$BP
    lead_snp <- lead_row$SNP
    lead_p   <- lead_row$P
    
    # Assign all SNPs within ±clump_bp of the lead SNP
    in_window <- (d$assigned == FALSE) &
      (d$BP >= (lead_bp - clump_bp)) &
      (d$BP <= (lead_bp + clump_bp))
    
    d[in_window, `:=`(
      assigned = TRUE,
      locus_id = locus_idx,
      lead_SNP = lead_snp,
      lead_BP  = lead_bp,
      lead_P   = lead_p
    )]
  }
  
  # One-row summary per locus
  d[, .(
    n_snps       = .N,
    minBP        = min(BP),
    maxBP        = max(BP),
    spanBP       = max(BP) - min(BP),
    lead_SNP     = lead_SNP[1],
    lead_BP      = lead_BP[1],
    lead_P       = lead_P[1],
    flank_kb     = clump_kb,
    window_start = pmax(lead_BP[1] - clump_bp, 1L),
    window_end   = lead_BP[1] + clump_bp
  ), by = .(locus_id)][order(lead_P)]
}

locus_tbl2 <- dt[, get_distance_loci(.SD), by = .(Element, CHR)]

# View the final recommended regions table
locus_tbl2[]

library(data.table)

fwrite(
  locus_tbl2,
  file = "/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_nine_main/locus_tbl2_distanceclump_flank500.csv"
)


#
fwrite(
  locus_tbl2[, .(Element, CHR, locus_id, lead_SNP, window_start, window_end)],
  file = "/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_nine_main/locus_tbl2_for_ld.tsv",
  sep = "\t"
)



#
library(data.table)
library(ggplot2)
library(ggrepel)

# -----------------------------
# File paths
# -----------------------------
gwas_file <- "/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_nine_main/Ca_no_outlier/Ca_no_outlier_MLMA.mlma"
ld_file   <- "/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_nine_main/Ca_chr4_locus1.ld"
gff_file  <- "/gpfs01/home/shxyt1/data/sticklebacks/ref/fresh/stickleback_v5_ensembl_genes.gff3.gz"

# -----------------------------
# Locus definition
# -----------------------------
lead_bp       <- 12857480L
chr_target    <- 4L
start_target  <- lead_bp - 250000L
end_target    <- lead_bp + 250000L
lead_target   <- "4:12857480"


# Label genes only within ±50 kb of the lead SNP
gene_label_start <- lead_bp - 50000L
gene_label_end   <- lead_bp + 50000L

# -----------------------------
# 1) Load FULL GWAS (GCTA-MLMA)
# -----------------------------
gwas_all <- fread(gwas_file, header = TRUE)
setnames(gwas_all, tolower(names(gwas_all)))  # chr, snp, bp, p, ...

# Keep SNPs in the plotting window
gwas_sub <- gwas_all[
  chr == chr_target & bp >= start_target & bp <= end_target
]

# Add plotting columns (ensure these always exist)
gwas_sub[, logp := -log10(p)]
gwas_sub[, is_lead := (snp == lead_target)]

# -----------------------------
# 2) Load LD (PLINK) and merge into GWAS
# -----------------------------
ld <- fread(ld_file)

# Map partner SNP -> r2 (PLINK output: SNP_A is lead, SNP_B is partner)
ld_map <- ld[SNP_A == lead_target, .(snp = SNP_B, r2 = R2)]
ld_map <- rbind(ld_map, data.table(snp = lead_target, r2 = 1), fill = TRUE)
ld_map <- unique(ld_map, by = "snp")

plot_dt <- merge(gwas_sub, ld_map, by = "snp", all.x = TRUE)

# LocusZoom-style r2 bins; NA -> grey
plot_dt[, r2_bin := fifelse(
  is.na(r2), "NA",
  fifelse(r2 >= 0.8, "0.8–1.0",
          fifelse(r2 >= 0.6, "0.6–0.8",
                  fifelse(r2 >= 0.4, "0.4–0.6",
                          fifelse(r2 >= 0.2, "0.2–0.4", "<0.2")
                  )
          )
  )
)]
plot_dt[, r2_bin := factor(r2_bin, levels = c("0.8–1.0","0.6–0.8","0.4–0.6","0.2–0.4","<0.2","NA"))]

# -----------------------------
# 3) Build genes_dt from GFF3 (no rtracklayer needed)
# -----------------------------
genes_dt <- fread(cmd = sprintf("zcat %s | grep -v '^#'", gff_file))
setnames(genes_dt, c("seqnames","source","type","start","end","score","strand","phase","attribute"))
genes_dt <- genes_dt[type == "gene"]

# Extract gene_id and gene_name from the GFF3 attribute column
genes_dt[, gene_id := sub(".*ID=([^;]+).*", "\\1", attribute)]
genes_dt[, gene_name := fifelse(
  grepl("Name=", attribute),
  sub(".*Name=([^;]+).*", "\\1", attribute),
  gene_id
)]

genes_dt <- genes_dt[, .(seqnames, start, end, strand, gene_id, gene_name)]

# Chromosome 4 in your GFF is chrIV
genes_chr4 <- genes_dt[seqnames == "chrIV"]
genes_win  <- genes_chr4[start <= end_target & end >= start_target]

# Genes to label: overlap ±50 kb around lead
genes_lab <- genes_chr4[start <= gene_label_end & end >= gene_label_start]
genes_lab[, label := fifelse(!is.na(gene_name) & gene_name != "", gene_name, gene_id)]

# Attach labels to genes in the plotting window (others keep NA label)
genes_win[, label := NA_character_]
if (nrow(genes_lab) > 0) {
  genes_win[gene_id %in% genes_lab$gene_id, label := genes_lab[label != "", label]]
}

# -----------------------------
# 4) Assign non-overlapping gene tracks (greedy)
# -----------------------------
setorder(genes_win, start, end)
genes_win[, track := 1L]

if (nrow(genes_win) > 1) {
  last_end <- integer(0)
  for (i in seq_len(nrow(genes_win))) {
    placed <- FALSE
    for (t in seq_len(max(genes_win$track))) {
      if (length(last_end) < t) last_end[t] <- -Inf
      if (genes_win$start[i] > last_end[t]) {
        genes_win$track[i] <- t
        last_end[t] <- genes_win$end[i]
        placed <- TRUE
        break
      }
    }
    if (!placed) {
      genes_win$track[i] <- max(genes_win$track) + 1L
      last_end[length(last_end) + 1L] <- genes_win$end[i]
    }
  }
}

# Limit gene rows shown (paper-friendly)
max_gene_rows <- 7L
genes_win <- genes_win[track <= max_gene_rows]
genes_win[, y := -track]

# -----------------------------
# 5) Build the plot (p)
# -----------------------------
ymax <- max(plot_dt$logp, na.rm = TRUE)
ymin <- -max_gene_rows - 0.8

# Assign an independent y to each gene's label (by track to offset layers)
setorder(genes_win, start)

n_tracks <- 3
genes_win[, track := ((.I - 1) %% n_tracks) + 1]

library(data.table)
library(ggplot2)

## =============================
## 1) Greedy packing helper
## =============================
pack_intervals <- function(dt, start_col, end_col, out_col = "track") {
  dt <- copy(dt)
  setorderv(dt, c(start_col, end_col))
  
  last_end <- numeric(0)
  track_id <- integer(nrow(dt))
  
  s <- dt[[start_col]]
  e <- dt[[end_col]]
  
  for (i in seq_len(nrow(dt))) {
    placed <- FALSE
    if (length(last_end) > 0) {
      for (t in seq_along(last_end)) {
        if (s[i] > last_end[t]) {
          track_id[i] <- t
          last_end[t] <- e[i]
          placed <- TRUE
          break
        }
      }
    }
    if (!placed) {
      last_end <- c(last_end, e[i])
      track_id[i] <- length(last_end)
    }
  }
  
  dt[, (out_col) := track_id]
  dt
}

## =============================
## 2) Prepare genes_win
## =============================
genes_tmp <- copy(genes_win)

# Label
genes_tmp[, label := fifelse(!is.na(gene_name) & gene_name != "", gene_name, gene_id)]

# Make very short genes visible
genes_tmp[, w := end - start]
genes_tmp[w < 20, `:=`(start = start - 10L, end = end + 10L)]
genes_tmp[, w := NULL]

## =============================
## 3) Body tracks (non-overlapping by gene interval)
## =============================
genes_tmp <- pack_intervals(genes_tmp, "start", "end", out_col = "body_track")

## =============================
## 4) Label tracks (non-overlapping by label width)
## =============================
region_bp <- end_target - start_target
bp_per_char <- region_bp / 80   # smaller denominator -> wider labels -> more tracks
pad_bp <- bp_per_char * 2

genes_tmp[, x_mid := (start + end) / 2]
genes_tmp[, lab_halfw := pmax(nchar(label) * bp_per_char / 2, bp_per_char * 2)]
genes_tmp[, lab_start := x_mid - lab_halfw - pad_bp]
genes_tmp[, lab_end   := x_mid + lab_halfw + pad_bp]

genes_tmp <- pack_intervals(genes_tmp, "lab_start", "lab_end", out_col = "lab_track")

# This is the final table for plotting
genes_lab <- genes_tmp

## =============================
## 5) Y positions (fixed near 0)
## =============================
ymax <- max(plot_dt$logp, na.rm = TRUE)

body_base_y <- -1.2
body_step   <- 0.75
genes_lab[, y_body := body_base_y - body_step * (body_track - 1L)]

lab_step   <- 0.55
lab_offset <- 0.22
lab_base_y <- min(genes_lab$y_body) - 0.35
genes_lab[, y_text := lab_base_y - lab_step * (lab_track - 1L) - lab_offset]

# Connector to show mapping
genes_lab[, cx := x_mid]
genes_lab[, cy := y_text]
genes_lab[, cyend := y_body]

ymin_plot <- min(genes_lab$y_text, na.rm = TRUE) - 0.8

## =============================
## 6) Plot (no ggrepel, because labels already non-overlapping)
## =============================
p <- ggplot() +
  
  geom_point(
    data = plot_dt[is_lead == FALSE],
    aes(x = bp, y = logp, color = r2_bin),
    size = 0.7, alpha = 0.45
  ) +
  
  geom_point(
    data = plot_dt[is_lead == TRUE],
    aes(x = bp, y = logp),
    size = 3.2, shape = 23, fill = "gold", color = "black", stroke = 0.4
  ) +
  
  geom_text(
    data = plot_dt[is_lead == TRUE],
    aes(x = bp, y = logp, label = snp),
    vjust = -0.8, size = 3
  ) +
  
  # Gene body
  geom_segment(
    data = genes_lab,
    aes(x = start, xend = end, y = y_body, yend = y_body),
    linewidth = 0.6, colour = "black"
  ) +
  
  # Connector
  geom_segment(
    data = genes_lab,
    aes(x = cx, xend = cx, y = cy, yend = cyend),
    linewidth = 0.25, colour = "grey40", alpha = 0.35
  ) +
  
  # Labels
  geom_text(
    data = genes_lab,
    aes(x = x_mid, y = y_text, label = label),
    size = 2.8
  ) +
  
  scale_color_manual(
    values = c(
      "0.8–1.0" = "#d73027",
      "0.6–0.8" = "#fc8d59",
      "0.4–0.6" = "#fee090",
      "0.2–0.4" = "#91bfdb",
      "<0.2"    = "#4575b4",
      "NA"      = "grey70"
    ),
    name = expression(LD~(r^2))
  ) +
  
  coord_cartesian(
    xlim = c(start_target, end_target),
    ylim = c(ymin_plot, ymax + 1)
  ) +
  
  labs(
    title = "Ca | chr4 locus1 (±500 kb around lead SNP)",
    x = "chr4 position (bp)",
    y = expression(-log[10](P))
  ) +
  
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

p



# 1) LD snps set (partners + lead)
ld_snps <- unique(c(ld$SNP_A, ld$SNP_B))

# 2) GWAS snps set in the plot window
gwas_snps <- unique(gwas_sub$snp)

# 3) overlap / coverage
cat("GWAS SNPs in window:", length(gwas_snps), "\n")
cat("LD SNPs in .ld file:", length(ld_snps), "\n")
cat("Overlap (GWAS ∩ LD):", sum(gwas_snps %in% ld_snps), "\n")
cat("Coverage in GWAS window:", round(mean(gwas_snps %in% ld_snps) * 100, 2), "%\n")

# 4) after merge, how many r2 are NA?
cat("After merge: total points:", nrow(plot_dt), "\n")
cat("After merge: r2 available:", sum(!is.na(plot_dt$r2)), "\n")
cat("After merge: r2 missing:", sum(is.na(plot_dt$r2)), "\n")
cat("Missing rate:", round(mean(is.na(plot_dt$r2)) * 100, 2), "%\n")

