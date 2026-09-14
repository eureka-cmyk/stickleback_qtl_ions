library(ggplot2)
options(bitmapType = "cairo")

pheno_dir <- "/gpfs01/home/shxyt1/data/sticklebacks/pheno/traits"

#  RINT =================== 
files <- list.files(
  pheno_dir,
  pattern = "^pheno_[A-Za-z]+\\.txt$",
  full.names = TRUE
)

for (file in files) {

  message("Processing: ", file)

  # pheno_Na.txt → Na
  element <- gsub("pheno_|\\.txt", "", basename(file))

  ph <- read.table(file, header = FALSE)
  colnames(ph) <- c("FID", "IID", element)

  ## ---------- 原始分布图 ----------
  out_png <- file.path(pheno_dir, paste0(element, "_distribution.png"))

  png(out_png, width = 1600, height = 1200, res = 150)
  print(
    ggplot(ph, aes_string(x = element)) +
      geom_histogram(bins = 40, fill = "steelblue", color = "black") +
      theme_bw() +
      ggtitle(paste("Distribution of", element))
  )
  dev.off()

  ## ---------- RINT ----------
  trait <- ph[[element]]
  rint <- qnorm((rank(trait, na.last = "keep") - 0.5) / sum(!is.na(trait)))
  ph[[paste0(element, "_rint")]] <- rint

  # 输出 RINT 文件
  out_rint <- file.path(pheno_dir, paste0("pheno_", element, "_RINT.txt"))
  write.table(
    ph[, c("FID", "IID", paste0(element, "_rint"))],
    out_rint,
    row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t"
  )

  ## ---------- RINT 分布图 ----------
  out_png_rint <- file.path(pheno_dir, paste0(element, "_distribution_RINT.png"))

  png(out_png_rint, width = 1600, height = 1200, res = 150)
  print(
    ggplot(ph, aes(x = rint)) +
      geom_histogram(bins = 40, fill = "steelblue", color = "black") +
      theme_bw() +
      ggtitle(paste("Distribution of", element, "after RINT"))
  )
  dev.off()

  message("Finished ", element, "\n")
}
#==========================

#=============log/zscore=============================================
library(data.table)

# Paths
base <- "/gpfs01/home/shxyt1/data/sticklebacks/pheno/traits"
orig_dir <- file.path(base, "original")
rint_dir <- file.path(base, "rint")
out_dir  <- file.path(base, "summary")

dir.create(out_dir, showWarnings = FALSE)

# Find original phenotype files

files <- list.files(
  orig_dir,
  pattern = "^pheno_[A-Za-z]+\\.txt$",
  full.names = TRUE
)

results <- list()

ph <- fread(
  "/gpfs01/home/shxyt1/data/sticklebacks/pheno/traits/original/pheno_Na.txt",
  header = FALSE
)

colnames(ph) <- c("FID", "IID", "VALUE")

ph[ph$IID == "QTL22149", ]

ph_log <- fread(
  "/gpfs01/home/shxyt1/data/sticklebacks/pheno/traits/log/Na_log.txt",
  header = FALSE
)

colnames(ph_log) <- c("FID", "IID", "VALUE")

ph_log[ph_log$IID == "QTL22149", ]

summary(ph_log$VALUE)


# Loop through each phenotype

for (file in files) {
  
  element <- gsub("pheno_|\\.txt", "", basename(file))
  message("Processing ", element)
  
  # ----- Read original -----
  ph <- fread(file, header = FALSE)
  colnames(ph) <- c("FID", "IID", "VALUE")
  
  raw <- ph$VALUE
  
  # Skip if all NA
  if (all(is.na(raw))) next
  
  # ----- log transform -----
  # Add epsilon to avoid log(0)
  log_val <- log(raw + 1e-6)
  
  # ----- z-score -----
  z_val <- scale(raw)[,1]
  
  # ----- read RINT transform -----
  rint_file <- file.path(rint_dir, paste0("pheno_", element, "_RINT.txt"))
  if (file.exists(rint_file)) {
    pr <- fread(rint_file, header = FALSE)
    rint_val <- pr$V3
  } else {
    rint_val <- rep(NA, length(raw))
  }
  
  # ===========================
  # Shapiro-Wilk tests
  # ===========================
  safe_shapiro <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) < 3) return(NA)
    if (length(x) > 5000) x <- sample(x, 5000)  # Shapiro limit
    shapiro.test(x)$p.value
  }
  
  p_raw    <- safe_shapiro(raw)
  p_log    <- safe_shapiro(log_val)
  p_zscore <- safe_shapiro(z_val)
  p_rint   <- safe_shapiro(rint_val)
  
  # Determine best transform
  all_ps <- c(raw = p_raw, log = p_log, zscore = p_zscore, rint = p_rint)
  best <- names(which.max(all_ps))
  
  # Save row
  results[[element]] <- data.frame(
    element = element,
    raw_p = p_raw,
    log_p = p_log,
    zscore_p = p_zscore,
    rint_p = p_rint,
    best_transform = best
  )
}


# Save summary table
summary_table <- rbindlist(results)
outfile <- file.path(out_dir, "summary_normality.csv")

fwrite(summary_table, outfile)

cat("Summary table saved to:\n", outfile, "\n")

#=================================

#===========log===========
library(data.table)

# Paths
base_dir <- "/gpfs01/home/shxyt1/data/sticklebacks/pheno/traits"

orig_dir <- file.path(base_dir, "original")
out_dir  <- file.path(base_dir, "log_no_outlier")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# original phenotype files
files <- list.files(
  orig_dir,
  pattern = "^pheno_[A-Za-z]+\\.txt$",
  full.names = TRUE
)

cat("Found", length(files), "phenotype files\n")

# Loop over traits
for (file in files) {
  
  element <- gsub("pheno_|\\.txt", "", basename(file))
  message("Processing ", element)
  
  # ----- Read original -----
  ph <- fread(file, header = FALSE)
  colnames(ph) <- c("FID", "IID", "VALUE")
  
  # ----- Pure log transform -----
  # IMPORTANT: no offset, no filtering, no tricks
  ph$LOG <- log(ph$VALUE)
  
  # ===========================
  # Write log-transformed data
  # (overwrite existing files)
  # ===========================
  out_data <- file.path(out_dir, paste0("pheno_", element, "_log.txt"))
  
  fwrite(
    ph[, .(FID, IID, LOG)],
    out_data,
    sep = "\t",
    col.names = FALSE
  )
  
  # ===========================
  # Distribution plot
  # ===========================
  out_plot <- file.path(out_dir, paste0(element, "_log_distribution.png"))
  
  png(out_plot, width = 1600, height = 1200, res = 150)
  print(
    ggplot(ph, aes(x = LOG)) +
      geom_histogram(bins = 40, fill = "steelblue", color = "black") +
      theme_bw() +
      ggtitle(paste("Log-transformed distribution of", element))
  )
  dev.off()
  
  message("Finished ", element, "\n")
}

cat("All log transformations completed.\n")
#===============================================

# log-delete-outlier=========
library(data.table)

#paths 
infile <- "/gpfs01/home/shxyt1/data/sticklebacks/pheno/traits/log_no_outlier/all_log_merged.txt"
outdir <- "/gpfs01/home/shxyt1/data/sticklebacks/pheno/traits/log_no_outlier"

#elements to process
elements <- c("As","Mg","Co","Cd","Sr","Pb","K","Fe","Ca")

# read data
dt <- fread(infile)

#function: IQR-based outlier detection 
is_outlier_iqr <- function(x, k = 3) {
  q1 <- quantile(x, 0.25, na.rm = TRUE)
  q3 <- quantile(x, 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  (x < (q1 - k * iqr)) | (x > (q3 + k * iqr))
}

#loop over elements 
for (el in elements) {
  
  pheno_col <- paste0("pheno_", el, "_log")
  stopifnot(pheno_col %in% colnames(dt))
  
  sub <- dt[, .(FID, IID, PHENO = get(pheno_col))]
  sub <- sub[!is.na(PHENO)]
  
  # detect outliers
  outlier_flag <- is_outlier_iqr(sub$PHENO)
  outliers <- sub[outlier_flag]
  clean <- sub[!outlier_flag]
  
  # save outliers 
  fwrite(
    outliers,
    file = file.path(outdir, paste0("outliers_", el, "_log.txt")),
    sep = "\t",
    col.names = TRUE
  )
  
  # save cleaned phenotype
  fwrite(
    clean,
    file = file.path(outdir, paste0("pheno_", el, "_no_outlier.txt")),
    sep = "\t",
    col.names = FALSE
  )
  
  cat(el, ": removed", nrow(outliers), "outliers; kept", nrow(clean), "samples\n")
}

library(data.table)

library(ggplot2)

indir  <- "/gpfs01/home/shxyt1/data/sticklebacks/pheno/traits/log_no_outlier"
outdir <- indir

elements <- c("As","Mg","Co","Cd","Sr","Pb","K","Fe","Ca")

for (el in elements) {
  
  infile <- file.path(indir, paste0("pheno_", el, "_no_outlier.txt"))
  dt <- fread(infile, col.names = c("FID","IID","PHENO"))
  
  p <- ggplot(dt, aes(x = PHENO)) +
    geom_histogram(aes(y = ..density..),
                   bins = 40,
                   fill = "#69b3a2",
                   color = "black",
                   alpha = 0.7) +
    geom_density(color = "red", linewidth = 1) +
    labs(
      title = paste0(el, " (log, no outlier)"),
      x = "Phenotype value",
      y = "Density"
    ) +
    theme_bw(base_size = 14)
  
  ggsave(
    filename = file.path(outdir, paste0(el, "_log_no_outlier_distribution.png")),
    plot = p,
    width = 6,
    height = 5,
    dpi = 300
  )
}

#PDF
indir  <- "/gpfs01/home/shxyt1/data/sticklebacks/pheno/traits/log_no_outlier"
outpdf <- file.path(indir, "log_no_outlier_distributions_9in1.pdf")

elements <- c("As","Mg","Co","Cd","Sr","Pb","K","Fe","Ca")

plots <- list()

for (el in elements) {
  
  infile <- file.path(indir, paste0("pheno_", el, "_no_outlier.txt"))
  dt <- fread(infile, col.names = c("FID","IID","PHENO"))
  
  p <- ggplot(dt, aes(x = PHENO)) +
    geom_histogram(aes(y = ..density..),
                   bins = 40,
                   fill = "#69b3a2",
                   color = "black",
                   alpha = 0.7) +
    geom_density(color = "red", linewidth = 0.8) +
    labs(
      title = el,
      x = NULL,
      y = "Density"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  
  plots[[el]] <- p
}

#combine 9 plots on one page 

combined <- wrap_plots(plots, ncol = 3)

# save to pdf 
ggsave(
  filename = outpdf,
  plot = combined,
  width = 12,
  height = 12
)

#====================================

#log-sig-nine-ele================
library(data.table)
library(ggplot2)
library(patchwork)

indir <- "/gpfs01/home/shxyt1/data/sticklebacks/pheno/traits/log_no_outlier"
outpdf <- file.path(indir, "significant_elements_distributions_2_col.pdf")

# ---------- element -> file mapping ----------
pheno_files <- list(
  Ca = "pheno_Ca_no_outlier.txt",
  Co = "pheno_Co_no_outlier.txt",
  Fe = "pheno_Fe_no_outlier.txt",
  Sr = "pheno_Sr_no_outlier.txt",
  P  = "pheno_P_log.txt",
  V  = "pheno_V_log.txt"
)

plots <- list()

for (el in names(pheno_files)) {
  
  infile <- file.path(indir, pheno_files[[el]])
  
  # no_outlier files: no header; log files: have header
  if (grepl("no_outlier", infile)) {
    dt <- fread(infile, col.names = c("FID","IID","PHENO"))
  } else {
    dt <- fread(infile)
    dt <- dt[, .(FID, IID, PHENO = get(paste0("pheno_", el, "_log")))]
  }
  
  dt <- dt[!is.na(PHENO)]
  
  p <- ggplot(dt, aes(x = PHENO)) +
    geom_histogram(aes(y = ..density..),
                   bins = 40,
                   fill = "#6BAED6",
                   color = "black",
                   alpha = 0.7) +
    geom_density(color = "red", linewidth = 0.8) +
    labs(
      title = el,
      x = NULL,
      y = "Density"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  
  plots[[el]] <- p
}

# ---------- combine 6 plots ----------
combined <- wrap_plots(plots, ncol = 2)

# ---------- save ----------
ggsave(
  filename = outpdf,
  plot = combined,
  width = 12,
  height = 12
)

#================================



#=========================box-cox===========
library(data.table)
library(ggplot2)
library(MASS)

# Paths
base_dir <- "/gpfs01/home/shxyt1/data/sticklebacks/pheno/traits"
orig_dir <- file.path(base_dir, "original")
out_dir  <- file.path(base_dir, "boxcox")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# original phenotype files
files <- list.files(
  orig_dir,
  pattern = "^pheno_[A-Za-z]+\\.txt$",
  full.names = TRUE
)

cat("Found", length(files), "phenotype files\n")

# Box-Cox transform function
boxcox_transform <- function(x, lambda) {
  if (abs(lambda) < 1e-8) return(log(x))
  (x^lambda - 1) / lambda
}

for (file in files) {
  
  element <- gsub("pheno_|\\.txt", "", basename(file))
  message("Processing ", element)
  
  # ----- Read original -----
  ph <- fread(file, header = FALSE)
  setnames(ph, c("FID", "IID", "VALUE"))
  
  # Ensure numeric
  ph[, VALUE := as.numeric(VALUE)]
  
  # Handle non-positive values by shifting (minimal shift to make min = small positive)
  minv <- min(ph$VALUE, na.rm = TRUE)
  shift <- 0
  if (minv <= 0) {
    shift <- abs(minv) + 1e-6
    message("  NOTE: non-positive values detected. Applying shift = ", shift)
  }
  
  x_pos <- ph$VALUE + shift
  
  # ----- Estimate lambda by MLE (intercept-only model) -----
  # boxcox() expects a fitted model; this finds lambda maximizing log-likelihood.
  fit <- lm(x_pos ~ 1)
  bc  <- boxcox(fit, lambda = seq(-2, 2, by = 0.01), plotit = FALSE)
  lambda_hat <- bc$x[which.max(bc$y)]
  
  # ----- Apply Box-Cox -----
  ph[, BOXCOX := boxcox_transform(x_pos, lambda_hat)]
  
  # ===========================
  # Write transformed data
  # ===========================
  out_data <- file.path(out_dir, paste0("pheno_", element, "_boxcox.txt"))
  
  fwrite(
    ph[, .(FID, IID, BOXCOX)],
    out_data,
    sep = "\t",
    col.names = FALSE
  )
  
  # ===========================
  # Distribution plot
  # ===========================
  out_plot <- file.path(out_dir, paste0(element, "_boxcox_distribution.png"))
  
  png(out_plot, width = 1600, height = 1200, res = 150)
  print(
    ggplot(ph, aes(x = BOXCOX)) +
      geom_histogram(bins = 40, fill = "steelblue", color = "black") +
      theme_bw() +
      ggtitle(paste0(
        "Box-Cox distribution of ", element,
        " (lambda=", round(lambda_hat, 3),
        if (shift > 0) paste0(", shift=", signif(shift, 3)) else "",
        ")"
      ))
  )
  dev.off()
  
  message("  lambda_hat = ", round(lambda_hat, 4))
  message("Finished ", element, "\n")
}

cat("All Box-Cox transformations completed.\n")






#=======combine data======
library(data.table)

boxcox_dir <- "/gpfs01/home/shxyt1/data/sticklebacks/pheno/traits/boxcox"
out_file <- file.path(boxcox_dir, "all_boxcox_merged.txt")

# list all *_boxcox.txt files
files <- list.files(
  boxcox_dir,
  pattern = "_boxcox\\.txt$",
  full.names = TRUE
)

cat("Found", length(files), "boxcox phenotype files\n")

# read first file as base
dt <- fread(files[1], header = FALSE)
setnames(dt, c("FID", "IID", gsub(".txt", "", basename(files[1]))))

# merge remaining files
for (f in files[-1]) {
  tmp <- fread(f, header = FALSE)
  colname <- gsub(".txt", "", basename(f))
  setnames(tmp, c("FID", "IID", colname))
  
  dt <- merge(dt, tmp, by = c("FID", "IID"), all = TRUE)
}

# write merged file (NO column names if you want GWAS style)
fwrite(
  dt,
  out_file,
  sep = "\t",
  col.names = TRUE
)

cat("Merged boxcox phenotype file written to:\n", out_file, "\n")
#======================

#============zscore==============
library(data.table)
library(ggplot2)

# Paths
base_dir <- "/gpfs01/home/shxyt1/data/sticklebacks/pheno/traits"

orig_dir <- file.path(base_dir, "original")
out_dir  <- file.path(base_dir, "zscore_no_outlier")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# original phenotype files
files <- list.files(
  orig_dir,
  pattern = "^pheno_[A-Za-z]+\\.txt$",
  full.names = TRUE
)

cat("Found", length(files), "phenotype files\n")

# Loop over traits
for (file in files) {
  
  element <- gsub("pheno_|\\.txt", "", basename(file))
  message("Processing ", element)
  
  # ----- Read original -----
  ph <- fread(file, header = FALSE)
  colnames(ph) <- c("FID", "IID", "VALUE")
  
  # ----- Z-score transform -----
  # (no filtering, NA will propagate naturally)
  mu <- mean(ph$VALUE, na.rm = TRUE)
  sdv <- sd(ph$VALUE, na.rm = TRUE)
  
  ph$ZSCORE <- (ph$VALUE - mu) / sdv
  
  # ===========================
  # Write z-score data
  # ===========================
  out_data <- file.path(out_dir, paste0("pheno_", element, "_zscore.txt"))
  
  fwrite(
    ph[, .(FID, IID, ZSCORE)],
    out_data,
    sep = "\t",
    col.names = FALSE
  )
  
  # ===========================
  # Distribution plot
  # ===========================
  out_plot <- file.path(out_dir, paste0(element, "_zscore_distribution.png"))
  
  png(out_plot, width = 1600, height = 1200, res = 150)
  print(
    ggplot(ph, aes(x = ZSCORE)) +
      geom_histogram(bins = 40, fill = "steelblue", color = "black") +
      theme_bw() +
      ggtitle(paste("Z-score distribution of", element))
  )
  dev.off()
  
  message("Finished ", element, "\n")
}

cat("All z-score transformations completed.\n")




#==============================

#==============combine distribution plot=====================
library(grid)
library(png)

# Paths
base_dir <- "/gpfs01/home/shxyt1/data/sticklebacks/pheno/traits/boxcox"
out_pdf  <- file.path(base_dir, "boxcox_distribution.pdf")

# Collect PNG files
png_files <- list.files(
  base_dir,
  pattern = "_boxcox_distribution\\.png$",
  full.names = TRUE
)

png_files <- sort(png_files)

cat("Found", length(png_files), "plots\n")

if (length(png_files) == 0) {
  stop("No PNG files found. Check directory or filename pattern.")
}

# PDF device
pdf(out_pdf, width = 16, height = 16)  # square page fits 4x4 nicely

# Plot pages
plots_per_page <- 16
n_pages <- ceiling(length(png_files) / plots_per_page)

idx <- 1

for (p in seq_len(n_pages)) {
  
  grid.newpage()
  
  pushViewport(
    viewport(
      layout = grid.layout(
        nrow = 4,
        ncol = 4
      )
    )
  )
  
  for (i in 1:plots_per_page) {
    
    if (idx > length(png_files)) break
    
    img <- readPNG(png_files[idx])
    
    row <- ((i - 1) %/% 4) + 1
    col <- ((i - 1) %% 4) + 1
    
    pushViewport(viewport(
      layout.pos.row = row,
      layout.pos.col = col
    ))
    
    grid.raster(img, interpolate = TRUE)
    
    popViewport()
    
    idx <- idx + 1
  }
  
  popViewport()
}

dev.off()

cat("PDF created:\n", out_pdf, "\n")


##========plot=============
#   Load libraries
library(data.table)
library(CMplot)
library(qqman)

# ===========================
# GWAS base directory
# ===========================
base_dir <- "/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_no_outlier"

dirs <- list.dirs(
  base_dir,
  recursive = FALSE,
  full.names = TRUE
)

cat("Found GWAS directories:\n")
print(dirs)

# ===========================
# Loop over all elements
# ===========================
for (d in dirs) {
  
  element <- basename(d)  # e.g. "Se_log"
  
  cat("\n==============================\n")
  cat("Processing:", element, "\n")
  cat("==============================\n")
  
  # MLMA file path: Se_log/Se_log_MLMA.mlma
  mlma_file <- file.path(d, paste0(element, "_MLMA.mlma"))
  
  if (!file.exists(mlma_file)) {
    cat("MLMA file not found — skip:\n", mlma_file, "\n")
    next
  }
  
  # Load MLMA results
  df <- fread(mlma_file, header = TRUE)
  
  # Prepare GWAS table for CMplot/qqman
  gwas <- data.frame(
    SNP = df$SNP,
    CHR = as.integer(df$Chr),
    BP  = df$bp,
    P   = df$p
  )
  
  gwas <- gwas[!is.na(gwas$CHR) & !is.na(gwas$BP) & !is.na(gwas$P), ]
  
  cat("Loaded SNP records:", nrow(gwas), "\n")
  
  # ===========================
  # Manhattan plot (CMplot)
  # Output goes into d
  # ===========================
  old_wd <- getwd()
  setwd(d)
  
  CMplot(
    gwas,
    plot.type = "m",
    LOG10 = TRUE,
    cex = 0.4,
    threshold = c(5e-8, 1e-6),
    threshold.col = c("red", "blue"),
    chr.den.col = NULL,
    signal.col = "red",
    file = "png",
    file.output = TRUE,
    dpi = 300,
    # element already contains "_log", so don't add another "_log"
    file.name = paste0(element),
    width = 14,
    height = 6
  )
  
  # ===========================
  # QQ plot (qqman)
  # Output goes into d
  # ===========================
  out_qq <- file.path(d, paste0(element, "_QQ.png"))
  
  png(out_qq, width = 1600, height = 1600, res = 200)
  qq(gwas$P, main = paste("QQ Plot —", element, "GWAS"))
  dev.off()
  
  setwd(old_wd)
  
  cat("Finished plots for:", element, "\n")
}





###with top snp===================
library(data.table)
base_dir <- "/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_nine_main"

dirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
dirs <- dirs[grepl("_no_outlier$", dirs)]


out_file <- "/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_nine_main/All_log_nine_mainTopSNPs.txt"
sink(out_file)

sig_threshold <- 5e-8
  
sugg_threshold <- 1e-6

cat("Element\tSNP\tCHR\tBP\tP\tCategory\n")

for (d in dirs) {
  
  element <- basename(d)
  prefix  <- sub("_no_outlier", "", element)
  
  mlma_file <- file.path(d, paste0(prefix, "_no_outlier_MLMA.mlma"))
  if (!file.exists(mlma_file)) {
    cat(prefix, "\tNA\tNA\tNA\tNA\tmissing_MLMA\n", sep="")
    next
  }
  
  df <- fread(mlma_file, header = TRUE)
  
  gwas <- data.frame(
    SNP = df$SNP,
    CHR = as.integer(df$Chr),
    BP  = df$bp,
    P   = df$p
  )
  gwas <- gwas[!is.na(gwas$CHR) & !is.na(gwas$BP) & !is.na(gwas$P), ]
  
  # Significant SNPs
  sig_snps <- gwas[gwas$P < sig_threshold, ]
  
  if (nrow(sig_snps) > 0) {
    for (i in 1:nrow(sig_snps)) {
      cat(prefix, sig_snps$SNP[i], sig_snps$CHR[i], sig_snps$BP[i], sig_snps$P[i],
          "significant", sep="\t")
      cat("\n")
    }
  }
  
  # Suggestive SNPs
  sugg_snps <- gwas[gwas$P < sugg_threshold & gwas$P >= sig_threshold, ]
  
  if (nrow(sugg_snps) > 0) {
    for (i in 1:nrow(sugg_snps)) {
      cat(prefix, sugg_snps$SNP[i], sugg_snps$CHR[i], sugg_snps$BP[i], sugg_snps$P[i],
          "suggestive", sep="\t")
      cat("\n")
    }
  }
  
  # No hits at all
  if (nrow(sig_snps) == 0 & nrow(sugg_snps) == 0) {
    cat(prefix, "\tNA\tNA\tNA\tNA\tno_hit\n", sep="")
  }
}

sink()
cat("All top SNPs saved to:\n", out_file, "\n")
