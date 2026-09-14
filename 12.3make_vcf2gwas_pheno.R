

infile  <- "/gpfs01/home/shxyt1/data/sticklebacks/pheno/all_pheno_clean.txt"
outfile <- "/gpfs01/home/shxyt1/data/sticklebacks/pheno/traits/vcf2gwas/vcf2gwas_pheno_Na.tsv"

ph <- read.table(infile, header = TRUE, sep = "\t", quote = "", comment.char = "", stringsAsFactors = FALSE, check.names = FALSE)


if (!("IID" %in% colnames(ph))) stop("No column named IID found in phenotype file.")

# ---- IID + Na ----
keep_cols <- c("IID", "Na")
missing_cols <- setdiff(keep_cols, colnames(ph))
if (length(missing_cols) > 0) stop(paste("Missing columns:", paste(missing_cols, collapse = ", ")))

out <- ph[, keep_cols]

# sample+trait
colnames(out) <- c("sample", "Na")

# data clean
out$Na <- suppressWarnings(as.numeric(out$Na))
out <- out[is.finite(out$Na), , drop = FALSE]

# remove duplication
out <- out[!duplicated(out$sample), , drop = FALSE]

# output
write.table(out, file = outfile, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

cat("Wrote:", outfile, "\n")
cat("N samples kept:", nrow(out), "\n")
