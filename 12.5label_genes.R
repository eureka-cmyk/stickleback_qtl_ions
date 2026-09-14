# install.packages("GenomicRanges")
# install.packages("BiocManager")
library(BiocManager)
BiocManager::install("rtracklayer", ask = FALSE, update = FALSE)
library(readr)
library(dplyr)
library(GenomicRanges)
library(IRanges)



# --- Step 1A: Read GWAS top SNPs (force column types) ---
gwas_file <- "/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_EDA10/All_EDA10_TopSNPs.txt"

df <- read_tsv(
  gwas_file,
  col_types = cols(
    Element  = col_character(),
    SNP      = col_character(),  # keep as character (avoid time parsing)
    CHR      = col_character(),
    BP       = col_integer(),
    P        = col_double(),
    Category = col_character()
  ),
  na = c("", "NA", ".")
) %>%
  mutate(
    CHR = gsub("^chr", "", CHR),
    BP  = as.integer(BP)
  ) %>%
  filter(!is.na(CHR), !is.na(BP))

# --- Step 1B: Quick sanity checks ---
cat("GWAS rows:", nrow(df), "\n")
cat("Unique CHR (first 25):", paste(head(sort(unique(df$CHR)), 25), collapse = ", "), "\n")
cat("Any NA in SNP?", any(is.na(df$SNP)), "\n")
cat("Any NA in BP?", any(is.na(df$BP)), "\n")
print(head(df, 3))

# --- Step 1C: Build variant GRanges (1bp loci) ---
vars_gr <- GRanges(
  seqnames = df$CHR,
  ranges   = IRanges(start = df$BP, end = df$BP),
  idx      = seq_len(nrow(df))
)

# Show a small preview
vars_gr[1:3]


library(data.table)
library(stringr)
library(GenomicRanges)
library(IRanges)

anno_file <- "/gpfs01/home/shxyt1/data/sticklebacks/ref/fresh/stickleback_v5_ensembl_genes.gff3.gz"

# --- Read GFF3 (skip comment lines starting with ##) ---
gff <- fread(
  cmd = paste("zcat", shQuote(anno_file), "| grep -v '^##'"),
  sep = "\t",
  header = FALSE,
  quote = "",
  data.table = FALSE
)

colnames(gff) <- c("seqid","source","type","start","end","score","strand","phase","attributes")

# --- Keep only gene features ---
genes <- gff[gff$type == "gene", c("seqid","start","end","strand","attributes")]
cat("Gene rows:", nrow(genes), "\n")

# --- Extract gene symbol from attributes ---
# Try common Ensembl GFF3 keys in order: gene_name, Name, gene, ID
extract_attr <- function(x, key) {
  m <- str_match(x, paste0("(^|;)", key, "=([^;]+)"))
  m[,3]
}

genes$gene_symbol <- extract_attr(genes$attributes, "gene_name")
if (all(is.na(genes$gene_symbol))) genes$gene_symbol <- extract_attr(genes$attributes, "Name")
if (all(is.na(genes$gene_symbol))) genes$gene_symbol <- extract_attr(genes$attributes, "gene")
if (all(is.na(genes$gene_symbol))) genes$gene_symbol <- extract_attr(genes$attributes, "ID")

cat("Non-NA gene_symbol:", sum(!is.na(genes$gene_symbol)), "\n")
print(head(genes[, c("seqid","start","end","gene_symbol")], 3))

# --- Build GRanges for genes ---
genes_gr <- GRanges(
  seqnames = genes$seqid,
  ranges   = IRanges(start = as.integer(genes$start), end = as.integer(genes$end)),
  strand   = genes$strand,
  gene_symbol = genes$gene_symbol
)

genes_gr[1:3]


library(dplyr)
library(stringr)

library(GenomeInfoDb)

# --- Check seqlevels before renaming ---
cat("Seqlevels in genes_gr2 BEFORE rename (all):\n")
# print(sort(seqlevels(genes_gr2)))
print(sort(seqlevels(genes_gr)))

# Also check if chrXXII exists
cat("Has chrXXII?", "chrXXII" %in% seqlevels(genes_gr2), "\n")


# --- Step 3A: Build mapping between GWAS CHR (1..22) and annotation seqid (chrI..chrXXII) ---
library(GenomicRanges)

# Roman -> Arabic mapping (names must match existing seqlevels)
roman <- c("I","II","III","IV","V","VI","VII","VIII","IX","X",
           "XI","XII","XIII","XIV","XV","XVI","XVII","XVIII","XIX","XX","XXI","XXII")

map <- setNames(as.character(1:22), paste0("chr", roman))

# 1) Check current seqlevels in genes_gr
old_levels <- seqlevels(genes_gr)
cat("Old seqlevels (first 30):", paste(head(old_levels, 30), collapse = ", "), "\n")

# 2) Keep only those we can map, drop the rest (e.g., chrM, scaffolds)
keep_old <- intersect(old_levels, names(map))
genes_gr2 <- keepSeqlevels(genes_gr, keep_old, pruning.mode = "coarse")

# 3) Rename seqlevels to Arabic (1..22)
genes_gr2 <- renameSeqlevels(genes_gr2, map)

# 4) Verify
cat("New seqlevels:", paste(sort(seqlevels(genes_gr2)), collapse = ", "), "\n")
cat("Example seqnames:", paste(head(as.character(seqnames(genes_gr2)), 5), collapse = ", "), "\n")


# --- Step 3B: Ensure genes_gr uses the same seqname style (chrI..chrXXII) ---
library(GenomicRanges)
library(IRanges)
library(dplyr)

# --- Step 3B-1: Split variants into annotatable (CHR 1-21) and non-annotatable ---
has_gene_chr <- sort(unique(as.character(seqlevels(genes_gr2))))  # "1".."21"

df_a <- df %>% filter(CHR %in% has_gene_chr)
df_b <- df %>% filter(!CHR %in% has_gene_chr)  # e.g., CHR==22

cat("Variants on annotated chrs:", nrow(df_a), "\n")
cat("Variants on unannotated chrs:", nrow(df_b), "\n")
cat("Unannotated CHR values:", paste(sort(unique(df_b$CHR)), collapse = ", "), "\n")

vars_a_gr <- GRanges(
  seqnames = df_a$CHR,
  ranges   = IRanges(start = df_a$BP, end = df_a$BP),
  idx      = seq_len(nrow(df_a))
)

# --- Step 3B-2: intragenic first ---
hits <- findOverlaps(vars_a_gr, genes_gr2, ignore.strand = TRUE)

intra_list <- split(mcols(genes_gr2)$gene_symbol[subjectHits(hits)], queryHits(hits))
intra_gene <- rep(NA_character_, length(vars_a_gr))
intra_gene[as.integer(names(intra_list))] <- vapply(
  intra_list,
  function(x) paste(sort(unique(x)), collapse = ";"),
  character(1)
)

# --- Step 3B-3: nearest for non-intragenic ---
no_intra <- which(is.na(intra_gene))

nearest_gene <- rep(NA_character_, length(vars_a_gr))
nearest_dist <- rep(NA_integer_, length(vars_a_gr))

if (length(no_intra) > 0) {
  nn <- distanceToNearest(vars_a_gr[no_intra], genes_gr2, ignore.strand = TRUE)
  nearest_gene[no_intra] <- mcols(genes_gr2)$gene_symbol[subjectHits(nn)]
  nearest_dist[no_intra] <- mcols(nn)$distance
}

# --- Step 3B-4: Assemble annotated part ---
out_a <- df_a %>%
  mutate(
    gene_symbol = if_else(!is.na(intra_gene), intra_gene, nearest_gene),
    gene_assign = if_else(!is.na(intra_gene), "intragenic", "nearest"),
    distance_bp = if_else(!is.na(intra_gene), 0L, nearest_dist)
  )

# --- Step 3B-5: Assemble unannotated part (e.g. chr22) ---
out_b <- df_b %>%
  mutate(
    gene_symbol = NA_character_,
    gene_assign = "unannotated_no_gene_chr",
    distance_bp = NA_integer_
  )

# --- Step 3B-6: Combine back (keep original row order as much as possible) ---
out <- bind_rows(out_a, out_b)

# Quick QC summary
cat("Intragenic:", sum(out$gene_assign == "intragenic"), "\n")
cat("Nearest:", sum(out$gene_assign == "nearest"), "\n")
cat("Unannotated:", sum(out$gene_assign == "unannotated_no_gene_chr"), "\n")
cat("Any NA gene_symbol overall:", sum(is.na(out$gene_symbol)), "\n")

print(head(out %>% dplyr::select(Element, SNP, CHR, BP, gene_symbol, gene_assign, distance_bp), 10))

# Save
write_csv(out, "/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_EDA10/All_EDA10_TopSNPs_with_genes.csv")


# --- Step 3C: intragenic first (overlap) ---
library(GenomicRanges)
library(IRanges)
library(dplyr)
library(readr)

# --- Step 3C-1: intragenic first ---
hits <- findOverlaps(vars_a_gr, genes_gr2, ignore.strand = TRUE)

intra_list <- split(mcols(genes_gr2)$gene_symbol[subjectHits(hits)], queryHits(hits))
intra_gene <- rep(NA_character_, length(vars_a_gr))
intra_gene[as.integer(names(intra_list))] <- vapply(
  intra_list,
  function(x) paste(sort(unique(x)), collapse = ";"),
  character(1)
)

# --- Step 3C-2: nearest for non-intragenic ---
no_intra <- which(is.na(intra_gene))

nearest_gene <- rep(NA_character_, length(vars_a_gr))
nearest_dist <- rep(NA_integer_, length(vars_a_gr))

if (length(no_intra) > 0) {
  nn <- distanceToNearest(vars_a_gr[no_intra], genes_gr2, ignore.strand = TRUE)
  nearest_gene[no_intra] <- mcols(genes_gr2)$gene_symbol[subjectHits(nn)]
  nearest_dist[no_intra] <- mcols(nn)$distance
}

# --- Step 3C-3: assemble outputs ---
out_a <- df_a %>%
  mutate(
    gene_symbol = if_else(!is.na(intra_gene), intra_gene, nearest_gene),
    gene_assign = if_else(!is.na(intra_gene), "intragenic", "nearest"),
    distance_bp = if_else(!is.na(intra_gene), 0L, nearest_dist)
  )

out_b <- df_b %>%
  mutate(
    gene_symbol = NA_character_,
    gene_assign = "unannotated_no_gene_chr",
    distance_bp = NA_integer_
  )

out <- bind_rows(out_a, out_b)

# --- Step 3C-4: summary checks ---
cat("Intragenic:", sum(out$gene_assign == "intragenic"), "\n")
cat("Nearest:", sum(out$gene_assign == "nearest"), "\n")
cat("Unannotated:", sum(out$gene_assign == "unannotated_no_gene_chr"), "\n")
cat("Any NA gene_symbol overall:", sum(is.na(out$gene_symbol)), "\n")

print(head(out %>% dplyr::select(Element, SNP, CHR, BP, gene_symbol, gene_assign, distance_bp), 10))

# --- Step 3C-5: save ---
out_file <- "/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_no_outlier/All_log_TopSNPs_with_genes.csv"
write_csv(out, out_file)
cat("Saved to:", out_file, "\n")

