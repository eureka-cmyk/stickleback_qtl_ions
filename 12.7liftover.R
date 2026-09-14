# ============================================================
# LiftOver QTL coordinates from v1 to v5 using a chain file
# Input:  qtl_catlog.xlsx (v1 coords)
# Output: qtl_catlog_lift_v1_to_v5.tsv (+ optional xlsx)
# ============================================================

library(data.table)
library(readxl)
library(rtracklayer)
library(GenomicRanges)

# ---- paths ----
in_xlsx     <- "/gpfs01/home/shxyt1/data/sticklebacks/qtl_catlog/qtl_catlog.xlsx"
chain_file  <- "/gpfs01/home/shxyt1/app/chain/v1_to_v5.chain.txt"
outdir      <- "/gpfs01/home/shxyt1/data/sticklebacks/qtl_catlog_v5"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---- helper: convert chrI/chrII/... to chr1/chr2/... if needed ----
roman_to_int <- function(x) {
  x <- gsub("^chr", "", x)
  x <- toupper(x)
  map <- c(I=1, II=2, III=3, IV=4, V=5, VI=6, VII=7, VIII=8, IX=9, X=10,
           XI=11, XII=12, XIII=13, XIV=14, XV=15, XVI=16, XVII=17, XVIII=18,
           XIX=19, XX=20, XXI=21)
  as.integer(map[x])
}

normalize_chr <- function(chr) {
  chr <- as.character(chr)
  chr[chr == "" | is.na(chr)] <- NA
  
  is_num   <- grepl("^chr[0-9]+$", chr)
  is_roman <- grepl("^chr[IVX]+$", chr)
  
  chr2 <- chr
  chr2[is_roman] <- paste0("chr", roman_to_int(chr2[is_roman]))
  chr2[is_num]   <- chr[is_num]
  chr2
}

# ---- helper: liftover one (chr, pos) pair ----
# Returns: new_chr, new_pos (v5). If failed -> NA.
lift_pair <- function(dt, chr_col, pos_col, chain) {
  
  chr <- normalize_chr(dt[[chr_col]])
  pos <- suppressWarnings(as.integer(dt[[pos_col]]))
  
  ok <- !is.na(chr) & !is.na(pos) & pos > 0
  
  new_chr <- rep(NA_character_, nrow(dt))
  new_pos <- rep(NA_integer_,   nrow(dt))
  
  if (!any(ok)) {
    return(list(new_chr = new_chr, new_pos = new_pos))
  }
  
  gr <- GRanges(seqnames = chr[ok], ranges = IRanges(start = pos[ok], end = pos[ok]))
  lo <- liftOver(gr, chain)
  
  idx_ok <- which(ok)
  for (i in seq_along(lo)) {
    if (length(lo[[i]]) >= 1) {
      picked <- lo[[i]][1]  # take the first mapping if multiple
      new_chr[idx_ok[i]] <- as.character(seqnames(picked))
      new_pos[idx_ok[i]] <- start(picked)
    }
  }
  
  list(new_chr = new_chr, new_pos = new_pos)
}

# ---- helper: liftover an interval defined by start and end ----
lift_interval4 <- function(dt, chr_s, pos_s, chr_e, pos_e, chain, prefix_out) {
  
  a <- lift_pair(dt, chr_s, pos_s, chain)
  b <- lift_pair(dt, chr_e, pos_e, chain)
  
  dt[[paste0(prefix_out, "_chr_start_v5")]] <- a$new_chr
  dt[[paste0(prefix_out, "_pos_start_v5")]] <- a$new_pos
  dt[[paste0(prefix_out, "_chr_end_v5")]]   <- b$new_chr
  dt[[paste0(prefix_out, "_pos_end_v5")]]   <- b$new_pos
  
  invisible(dt)
}

# ---- helper: only liftover if all required columns exist ----
maybe_lift <- function(dt, chr_s, pos_s, chr_e, pos_e, prefix_out, chain) {
  need <- c(chr_s, pos_s, chr_e, pos_e)
  if (all(need %in% names(dt))) {
    message("Lifting: ", prefix_out, " (", paste(need, collapse = ", "), ")")
    lift_interval4(dt, chr_s, pos_s, chr_e, pos_e, chain, prefix_out)
  } else {
    message("Skip (missing cols): ", prefix_out)
  }
}

# ---- load data ----
dt <- as.data.table(read_excel(in_xlsx, sheet = 1))

# ---- load chain ----
ch <- import.chain(chain_file)

# ---- liftover key blocks ----
maybe_lift(dt, "BS1_chr_start",     "BS1_pos_start",     "BS1_chr_end",     "BS1_pos_end",     "BS1_marker", ch)
maybe_lift(dt, "BS1_chr_regstart",  "BS1_pos_regstart",  "BS1_chr_regend",  "BS1_pos_regend",  "BS1_region", ch)

maybe_lift(dt, "ROE_chr_start",     "ROE_pos_start",     "ROE_chr_end",     "ROE_pos_end",     "ROE_marker", ch)
maybe_lift(dt, "ROE_chr_regstart",  "ROE_pos_regstart",  "ROE_chr_regend",  "ROE_pos_regend",  "ROE_region", ch)

maybe_lift(dt, "GLA_chr_start",     "GLA_pos_start",     "GLA_chr_end",     "GLA_pos_end",     "GLA_marker", ch)
maybe_lift(dt, "GLA_chr_regstart",  "pos_regstart_GLA",  "GLA_chr_regend",  "GLA_pos_regend",  "GLA_region", ch)

# ---- save ----
out_tsv  <- file.path(outdir, "qtl_catlog_lift_v1_to_v5.tsv")
fwrite(dt, out_tsv, sep = "\t", quote = FALSE, na = "NA")
cat("Saved:", out_tsv, "\n")

# Optional: also save xlsx
if (requireNamespace("openxlsx", quietly = TRUE)) {
  out_xlsx <- file.path(outdir, "qtl_catlog_lift_v1_to_v5.xlsx")
  openxlsx::write.xlsx(as.data.frame(dt), out_xlsx, overwrite = TRUE)
  cat("Saved:", out_xlsx, "\n")
}




#
suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(GenomicRanges)
})

#===================== 0) Paths =====================
GWAS_CSV <- "/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_nine_main/Six_log_main_TopSNPs_with_genes.csv"
QTL_XLSX <- "/gpfs01/home/shxyt1/data/sticklebacks/qtl_catlog/qtl_catlog.v1_to_v5.lifted.xlsx"

OUT_TSV  <- "/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_nine_main/Six_log_main_TopSNPs_with_QTLcatalog.tsv"

#===================== 1) Load GWAS =====================
gwas <- fread(GWAS_CSV)

# Keep significant only (adjust if you want all rows)
if ("Category" %in% names(gwas)) {
  gwas <- gwas[Category == "significant"]
}

# Ensure types
gwas[, CHR := as.integer(CHR)]
gwas[, BP  := as.integer(BP)]

# Map 1-21 -> chrI..chrXXI (to match QTL chr naming)
rom <- c("chrI","chrII","chrIII","chrIV","chrV","chrVI","chrVII","chrVIII","chrIX","chrX",
         "chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI","chrXVII","chrXVIII","chrXIX","chrXX","chrXXI")
gwas[, chr := ifelse(!is.na(CHR) & CHR>=1 & CHR<=21, rom[CHR], NA_character_)]

# Drop rows without chr/BP
gwas <- gwas[!is.na(chr) & !is.na(BP) & BP > 0]
gwas[, gwas_id := .I]  # stable id for join back

#===================== 2) Load QTL catalog (v5-lifted) =====================
# IMPORTANT: your original xlsx has a "big header row", but the lifted.xlsx we wrote has a clean header in row 1.
# If for any reason you still have the big-header format, change col_names=TRUE & try skip=1.
qtl <- as.data.table(read_excel(QTL_XLSX, sheet = 1, col_names = TRUE))

# Required metadata columns you want to attach
need_meta <- c("ref","populations","ecotypes","species","qtl_trait","qtl_cat")
miss_meta <- setdiff(need_meta, names(qtl))
if (length(miss_meta) > 0) {
  stop("QTL table missing metadata columns: ", paste(miss_meta, collapse=", "))
}

#===================== 3) Build union QTL interval in v5 (scheme A) =====================
# We use REGION coordinates (recommended). For each QTL row, union across BS1/ROE/GLA region in v5.

# Helper: pick first non-NA among multiple columns
first_non_na <- function(...) {
  xs <- list(...)
  out <- rep(NA_character_, length(xs[[1]]))
  for (x in xs) {
    idx <- is.na(out) & !is.na(x) & x != ""
    out[idx] <- as.character(x[idx])
  }
  out
}

# Ensure numeric
num_cols <- c("BS1_region_pos_start_v5","BS1_region_pos_end_v5",
              "ROE_region_pos_start_v5","ROE_region_pos_end_v5",
              "GLA_region_pos_start_v5","GLA_region_pos_end_v5")
for (cc in intersect(num_cols, names(qtl))) {
  qtl[, (cc) := as.integer(get(cc))]
}

# Chrom columns (region chr start in v5)
chr_cols <- c("BS1_region_chr_start_v5","ROE_region_chr_start_v5","GLA_region_chr_start_v5")
for (cc in intersect(chr_cols, names(qtl))) {
  qtl[, (cc) := as.character(get(cc))]
}

# If any of these v5 region columns are missing, stop (because scheme A relies on them)
need_cols <- c(chr_cols, num_cols)
miss_cols <- setdiff(need_cols, names(qtl))
if (length(miss_cols) > 0) {
  stop("QTL lifted table missing expected v5 region columns: ", paste(miss_cols, collapse=", "),
       "\nTip: ensure you generated *_region_chr_*_v5 and *_region_pos_*_v5 in the lift+backfill step.")
}

# Union chr: first available among BS1/ROE/GLA region chr
qtl[, union_chr := first_non_na(BS1_region_chr_start_v5, ROE_region_chr_start_v5, GLA_region_chr_start_v5)]

# Union start/end: min start across available; max end across available (row-wise)
qtl[, union_start := do.call(pmin, c(list(BS1_region_pos_start_v5, ROE_region_pos_start_v5, GLA_region_pos_start_v5), na.rm=TRUE))]
qtl[, union_end   := do.call(pmax, c(list(BS1_region_pos_end_v5,   ROE_region_pos_end_v5,   GLA_region_pos_end_v5),   na.rm=TRUE))]

# Clean: invalid rows
qtl <- qtl[!is.na(union_chr) & !is.na(union_start) & !is.na(union_end) & union_start > 0 & union_end >= union_start]

# Make a QTL id
qtl[, qtl_id := .I]

#===================== 4) Overlap: GWAS SNP in union QTL interval =====================
gr_gwas <- GRanges(seqnames = gwas$chr, ranges = IRanges(start = gwas$BP, end = gwas$BP), gwas_id = gwas$gwas_id)

gr_qtl  <- GRanges(seqnames = qtl$union_chr,
                   ranges = IRanges(start = qtl$union_start, end = qtl$union_end),
                   qtl_id = qtl$qtl_id)

hits <- findOverlaps(gr_gwas, gr_qtl, ignore.strand = TRUE)

# If no overlap at all, still write output with empty annotations
gwas[, previous_qtl := "No"]
gwas[, `:=`(
  qtl_ref = NA_character_,
  qtl_populations = NA_character_,
  qtl_ecotypes = NA_character_,
  qtl_species = NA_character_,
  qtl_trait = NA_character_,
  qtl_cat = NA_character_
)]

if (length(hits) > 0) {
  hit_dt <- data.table(
    gwas_id = mcols(gr_gwas)$gwas_id[queryHits(hits)],
    qtl_id  = mcols(gr_qtl)$qtl_id[subjectHits(hits)]
  )
  
  # Attach metadata
  hit_dt <- merge(hit_dt, qtl[, .(qtl_id, ref, populations, ecotypes, species, qtl_trait, qtl_cat)], by="qtl_id", all.x=TRUE)
  
  # Collapse multiple QTL matches per GWAS SNP
  agg <- hit_dt[, .(
    qtl_ref         = paste(unique(ref), collapse="; "),
    qtl_populations = paste(unique(populations), collapse="; "),
    qtl_ecotypes    = paste(unique(ecotypes), collapse="; "),
    qtl_species     = paste(unique(species), collapse="; "),
    qtl_trait       = paste(unique(qtl_trait), collapse="; "),
    qtl_cat         = paste(unique(qtl_cat), collapse="; ")
  ), by = gwas_id]
  
  # Join back
  gwas <- merge(gwas, agg, by="gwas_id", all.x=TRUE, sort=FALSE)
  
  # Set flag (robust: based on overlap ids, not column names after merge)
  gwas[gwas_id %in% agg$gwas_id, previous_qtl := "Yes"]
}

# Restore original order
setorder(gwas, gwas_id)

#===================== 5) Save =====================
fwrite(gwas, OUT_TSV, sep="\t", quote=FALSE, na="NA")
cat("Saved:", OUT_TSV, "\n")
cat("GWAS SNPs:", nrow(gwas), "\n")
cat("Overlapping (previous_qtl=Yes):", sum(gwas$previous_qtl=="Yes"), "\n")
# 
