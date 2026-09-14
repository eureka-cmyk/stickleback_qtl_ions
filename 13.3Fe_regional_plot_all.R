
library(dplyr)

# 1. Define Paths (Updated to your requested directory)
input_file  <- "/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_EDA10/All_EDA10_TopSNPs_with_genes.csv"
output_dir  <- "/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_EDA10"
output_file <- file.path(output_dir, "locus_tbl2_for_ld.tsv")

# 2. Set Parameters (You can adjust these)
# merge_dist: If two SNPs are within this distance, they are considered one locus
merge_dist <- 500000  
# window_size: Total buffer around the lead SNP
window_buffer <- 500000 

# 3. Read Data
raw_data <- read.csv(input_file, stringsAsFactors = FALSE)

# 4. Process Data
locus_table <- raw_data %>%
  filter(Category == "significant") %>%
  mutate(Element = sub("_.*", "", Trait)) %>%
  arrange(Element, CHR, BP) %>%
  group_by(Element, CHR) %>%
  # Logic to merge nearby SNPs into the same Locus ID
  mutate(
    dist_to_prev = BP - lag(BP, default = 0),
    is_new_locus = if_else(dist_to_prev > merge_dist, 1, 0),
    locus_id = cumsum(is_new_locus)
  ) %>%
  group_by(Element, CHR, locus_id) %>%
  # For each locus, pick the most significant SNP (lowest P-value) as the lead
  summarize(
    lead_SNP = SNP[which.min(P)],
    lead_BP = BP[which.min(P)],
    .groups = "drop"
  ) %>%
  # Define the final window based on the Lead SNP position
  mutate(
    window_start = pmax(1, lead_BP - window_buffer),
    window_end = lead_BP + window_buffer
  ) %>%
  select(Element, CHR, locus_id, lead_SNP, window_start, window_end)

# 5. Save Output
write.table(locus_table, 
            file = output_file, 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE)

cat("File saved to:", output_file, "\n")
cat("Total loci identified after merging:", nrow(locus_table), "\n")




# --- 1. Path Settings ---
# Input: Your original GWAS top SNPs file
input_file  <- "/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_EDA10/All_EDA10_TopSNPs_with_genes.csv"
output_dir  <- "/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_EDA10"

# Output: Saving to a NEW name to avoid overwriting your existing tsv
output_csv  <- file.path(output_dir, "locus_summary_detailed_stats.csv")

# --- 2. Parameter Settings ---
merge_dist <- 500000   # Distance to group SNPs into one locus (500kb)
flank_bp   <- 500000   # Flank for window_start/end (500kb)

# --- 3. Processing ---
raw_data <- read.csv(input_file, stringsAsFactors = FALSE)

detailed_locus_table <- raw_data %>%
  # Only significant SNPs
  filter(Category == "significant") %>%
  # Extract Element (e.g., Ca, Co, Fe)
  mutate(Element = sub("_.*", "", Trait)) %>%
  # Sort to ensure clustering works correctly
  arrange(Element, CHR, BP) %>%
  # Grouping Logic: identify clusters of SNPs
  group_by(Element, CHR) %>%
  mutate(
    is_new_locus = if_else(BP - lag(BP, default = first(BP)) > merge_dist, 1, 0),
    locus_id = cumsum(is_new_locus) + 1
  ) %>%
  # Calculate the 13 columns of interest
  group_by(Element, CHR, locus_id) %>%
  summarize(
    n_snps   = n(),
    minBP    = min(BP),
    maxBP    = max(BP),
    spanBP   = maxBP - minBP,
    lead_SNP = SNP[which.min(P)],
    lead_BP  = BP[which.min(P)],
    lead_P   = min(P),
    .groups  = "drop"
  ) %>%
  # Add window and flank info
  mutate(
    flank_kb     = flank_bp / 1000,
    window_start = pmax(1, lead_BP - flank_bp),
    window_end   = lead_BP + flank_bp
  )

# --- 4. Save Output ---
# This saves as a comma-separated CSV as you requested
write.csv(detailed_locus_table, 
          file = output_csv, 
          row.names = FALSE, 
          quote = FALSE)

cat("Done! The detailed CSV has been saved to:", output_csv, "\n")

## 0) Paths / parameters
## =============================
# Locus table (contains Element/CHR/locus_id/lead_SNP/lead_BP)
locus_file <- "/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_EDA10/locus_summary_detailed_stats.csv"

# GFF3 for gene tracks
gff_file <- "/gpfs01/home/shxyt1/data/sticklebacks/ref/fresh/stickleback_v5_ensembl_genes.gff3.gz"

# LD directory (LD files; Fe only will be used)
ld_dir <- "/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_EDA10/LD_all"

# Output directory (auto-create)
out_dir <- "/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_EDA10/regional plot"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Plot window settings: ±150 kb
plot_flank_bp <- 500000L   # ±150 kb
gene_label_flank_bp <- 500000L  # label genes within ±100 kb of lead (you can change)

# PNG size
png_width  <- 800
png_height <- 588
png_res    <- 96

## =============================
## 1) Trait -> GWAS file mapping
## =============================
get_gwas_path <- function(trait) {
  trait <- as.character(trait)
  if (trait %in% c("Ca","Co","Fe","Sr")) {
    return(sprintf(
      "/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_EDA10/%s_no_outlier/%s_no_outlier_MLMA.mlma",
      trait, trait
    ))
  }
  if (trait %in% c("P","V")) {
    return(sprintf(
      "/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_EDA10/%s_log/%s_log_MLMA.mlma",
      trait, trait
    ))
  }
  stop("Unknown trait: ", trait)
}

## =============================
## 2) chr number -> GFF seqname (chrI, chrII, ...)
## =============================
to_roman <- function(n) {
  romans <- c("I","II","III","IV","V","VI","VII","VIII","IX","X",
              "XI","XII","XIII","XIV","XV","XVI","XVII","XVIII","XIX","XX","XXI",
              "XXII","XXIII","XXIV","XXV")
  if (n < 1 || n > length(romans)) return(NA_character_)
  romans[n]
}
chr_to_gff <- function(chr_int) {
  r <- to_roman(as.integer(chr_int))
  if (is.na(r)) return(NA_character_)
  paste0("chr", r)
}

## =============================
## 3) Load gene annotation ONCE
## =============================
message("Loading GFF genes once...")
genes_dt <- fread(cmd = sprintf("zcat %s | grep -v '^#'", gff_file))
setnames(genes_dt, c("seqnames","source","type","start","end","score","strand","phase","attribute"))
genes_dt <- genes_dt[type == "gene"]

genes_dt[, gene_id := sub(".*ID=([^;]+).*", "\\1", attribute)]
genes_dt[, gene_name := fifelse(
  grepl("Name=", attribute),
  sub(".*Name=([^;]+).*", "\\1", attribute),
  gene_id
)]
genes_dt <- genes_dt[, .(seqnames, start=as.integer(start), end=as.integer(end),
                         strand, gene_id, gene_name)]
message("GFF loaded: ", nrow(genes_dt), " genes")

## =============================
## 4) Non-overlapping packing helper
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
## 5) One-locus plotting function
## =============================
plot_one_locus <- function(trait, chr_target, locus_id, lead_snp, lead_bp) {
  
  # GWAS path
  gwas_file <- get_gwas_path(trait)
  if (!file.exists(gwas_file)) {
    warning("GWAS file not found: ", gwas_file)
    return(FALSE)
  }
  
  # windows
  chr_target <- as.integer(chr_target)
  lead_bp <- as.integer(lead_bp)
  start_target <- max(1L, lead_bp - plot_flank_bp)
  end_target   <- lead_bp + plot_flank_bp
  
  gene_label_start <- max(1L, lead_bp - gene_label_flank_bp)
  gene_label_end   <- lead_bp + gene_label_flank_bp
  
  # 1) Load GWAS and subset
  gwas_all <- fread(gwas_file, header = TRUE)
  setnames(gwas_all, tolower(names(gwas_all)))  # expect chr, snp, bp, p
  
  if (!all(c("chr","snp","bp","p") %in% names(gwas_all))) {
    warning("GWAS columns missing in: ", gwas_file, " | got: ", paste(names(gwas_all), collapse=", "))
    return(FALSE)
  }
  
  gwas_sub <- gwas_all[chr == chr_target & bp >= start_target & bp <= end_target]
  if (nrow(gwas_sub) == 0) {
    warning(sprintf("%s chr%d locus%d: no SNPs in window [%d, %d].",
                    trait, chr_target, locus_id, start_target, end_target))
    return(FALSE)
  }
  
  gwas_sub[, logp := -log10(p)]
  gwas_sub[, is_lead := (snp == lead_snp)]
  
  # 2) Load LD and merge (if exists)
  ld_file <- file.path(ld_dir, sprintf("%s_chr%d_locus%d_lead%s.ld",
                                       trait, chr_target, locus_id, lead_snp))
  if (file.exists(ld_file)) {
    ld <- fread(ld_file)
    if (!all(c("SNP_A","SNP_B","R2") %in% names(ld))) {
      warning("LD file missing SNP_A/SNP_B/R2: ", ld_file)
      ld_map <- data.table(snp = character(), r2 = numeric())
    } else {
      ld_map <- ld[SNP_A == lead_snp, .(snp = SNP_B, r2 = as.numeric(R2))]
      ld_map <- rbind(ld_map, data.table(snp = lead_snp, r2 = 1), fill = TRUE)
      ld_map <- unique(ld_map, by = "snp")
    }
  } else {
    ld_map <- data.table(snp = character(), r2 = numeric())
  }
  
  plot_dt <- merge(gwas_sub, ld_map, by = "snp", all.x = TRUE)
  
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
  plot_dt[, r2_bin := factor(r2_bin,
                             levels = c("0.8–1.0","0.6–0.8","0.4–0.6","0.2–0.4","<0.2","NA"))]
  
  # 3) Genes
  gff_seq <- chr_to_gff(chr_target)
  genes_chr <- genes_dt[seqnames == gff_seq]
  genes_win <- genes_chr[start <= end_target & end >= start_target]
  
  if (nrow(genes_win) > 0) {
    genes_win <- copy(genes_win)
    
    genes_win[, label := fifelse(!is.na(gene_name) & gene_name != "", gene_name, gene_id)]
    genes_win[, w := end - start]
    genes_win[w < 20, `:=`(start = start - 10L, end = end + 10L)]
    genes_win[, w := NULL]
    
    # Label only those overlapping ±gene_label_flank_bp around lead
    genes_win[, label := fifelse(start <= gene_label_end & end >= gene_label_start, label, NA_character_)]
    
    # Tracks for gene bodies
    genes_win <- pack_intervals(genes_win, "start", "end", out_col = "body_track")
    
    # Tracks for labels based on label widths
    region_bp <- end_target - start_target
    bp_per_char <- region_bp / 80
    pad_bp <- bp_per_char * 2
    
    genes_win[, x_mid := (start + end) / 2]
    genes_win[, lab_halfw := pmax(nchar(fifelse(is.na(label), "", label)) * bp_per_char / 2, bp_per_char * 2)]
    genes_win[, lab_start := x_mid - lab_halfw - pad_bp]
    genes_win[, lab_end   := x_mid + lab_halfw + pad_bp]
    
    genes_lab <- genes_win[!is.na(label)]
    if (nrow(genes_lab) > 0) {
      genes_lab <- pack_intervals(genes_lab, "lab_start", "lab_end", out_col = "lab_track")
    } else {
      genes_lab <- genes_win[0]
      genes_lab[, lab_track := integer()]
    }
    
    # y positions
    body_base_y <- -1.2
    body_step   <- 0.75
    genes_win[, y_body := body_base_y - body_step * (body_track - 1L)]
    
    if (nrow(genes_lab) > 0) {
      lab_step   <- 0.55
      lab_offset <- 0.22
      lab_base_y <- min(genes_win$y_body) - 0.35
      
      genes_lab[, y_text := lab_base_y - lab_step * (lab_track - 1L) - lab_offset]
      
      # connector
      genes_lab[, cx := x_mid]
      genes_lab[, cy := y_text]
      y_map <- unique(genes_win[, .(gene_id, y_body)])
      genes_lab <- merge(genes_lab, y_map, by = "gene_id", all.x = TRUE)
      
      if (any(is.na(genes_lab$y_body))) {
        genes_lab[is.na(y_body), y_body := min(genes_win$y_body, na.rm = TRUE)]
      }
      
      genes_lab[, cyend := y_body]
      ymin_plot <- min(genes_lab$y_text, na.rm = TRUE) - 0.8
    } else {
      ymin_plot <- min(genes_win$y_body, na.rm = TRUE) - 0.8
    }
  } else {
    genes_win <- data.table()
    genes_lab <- data.table()
    ymin_plot <- -1.8
  }
  
  ymax <- max(plot_dt$logp, na.rm = TRUE)
  
  # 4) Build plot
  title_txt <- sprintf("%s|chr%dlocus%d|±75kb", trait, chr_target, locus_id)
  
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
    {if (nrow(genes_win) > 0) geom_segment(
      data = genes_win,
      aes(x = start, xend = end, y = y_body, yend = y_body),
      linewidth = 0.6, colour = "black"
    )} +
    {if (nrow(genes_lab) > 0) geom_segment(
      data = genes_lab,
      aes(x = cx, xend = cx, y = cy, yend = cyend),
      linewidth = 0.25, colour = "grey40", alpha = 0.35
    )} +
    {if (nrow(genes_lab) > 0) geom_text(
      data = genes_lab,
      aes(x = x_mid, y = y_text, label = label),
      size = 2.8
    )} +
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
      title = title_txt,
      x = sprintf("chr%d position (bp)", chr_target),
      y = expression(-log[10](P))
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "right",
      panel.grid.minor = element_blank()
    )
  
  # 5) Save PNG (include _150kb in name)
  safe_name <- sprintf("%s_chr%d_locus%d_150kb.png", trait, chr_target, locus_id)
  out_png <- file.path(out_dir, safe_name)
  
  png(out_png, width = png_width, height = png_height, res = png_res, type = "cairo-png")
  print(p)
  dev.off()
  
  message("Saved: ", out_png, if (file.exists(ld_file)) "" else " (LD missing -> grey)")
  TRUE
}

## =============================
## 6) Batch run (Fe only)
## =============================
locus_tbl <- fread(locus_file)

need_cols <- c("Element","CHR","locus_id","lead_SNP","lead_BP")
miss <- setdiff(need_cols, names(locus_tbl))
if (length(miss) > 0) stop("Missing columns in locus table: ", paste(miss, collapse=", "))

# Fe only
locus_tbl <- locus_tbl[Element == "Fe"]
if (nrow(locus_tbl) == 0) stop("No Fe rows found in locus table!")

locus_tbl[, CHR := as.integer(CHR)]
locus_tbl[, locus_id := as.integer(locus_id)]
locus_tbl[, lead_BP := as.integer(lead_BP)]
locus_tbl[, Element := as.character(Element)]
locus_tbl[, lead_SNP := as.character(lead_SNP)]

ok <- logical(nrow(locus_tbl))
for (i in seq_len(nrow(locus_tbl))) {
  tr   <- locus_tbl$Element[i]   # Fe
  ch   <- locus_tbl$CHR[i]
  li   <- locus_tbl$locus_id[i]
  lsnp <- locus_tbl$lead_SNP[i]
  lbp  <- locus_tbl$lead_BP[i]
  
  message(sprintf("Plotting %s chr%d locus%d lead %s (%d) ...", tr, ch, li, lsnp, lbp))
  ok[i] <- plot_one_locus(tr, ch, li, lsnp, lbp)
}

message("Done. Success: ", sum(ok), " / ", length(ok))

