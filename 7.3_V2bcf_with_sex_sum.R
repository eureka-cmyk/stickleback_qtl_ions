library(tidyverse)
library(patchwork)
library(viridis)

options(bitmapType = "cairo")

# ============================
# 1. Define paths and settings
# ============================

# New VCFtools QC results
stat_dir <- "/gpfs01/home/shxyt1/data/sticklebacks/VCF/stats_sex_corrected"

# Sample metadata
meta_file <- "/gpfs01/home/shxyt1/data/sticklebacks/pheno/qtl_all_phenotype_data.csv"

# Sample ID column in metadata
sample_col <- "fish_id"

# Grouping variable for plots and heterozygosity comparison
group_var <- "sex"

# QC cutoffs
missing_cutoff <- 0.40
excess_het_cutoff <- 0.015

# Sex chromosomes excluded from routine sample QC
excluded_chr <- c(
  "NC_053230.1",  # chromosome 19 / X
  "NC_053233.1"   # Y chromosome
)

# ============================
# 2. List chromosome prefixes
# ============================

het_files_all <- list.files(
  stat_dir,
  pattern = "\\.het$",
  full.names = TRUE
)

if (length(het_files_all) == 0) {
  stop("No .het files found in: ", stat_dir)
}

# Remove .het suffix
prefixes_all <- sub("\\.het$", "", het_files_all)

# Exclude X/chr19 and Y chromosome
prefixes <- prefixes_all[
  !str_detect(
    basename(prefixes_all),
    str_c(excluded_chr, collapse = "|")
  )
]

cat("All chromosome files found:", length(prefixes_all), "\n")
cat("Autosomal chromosome files retained:", length(prefixes), "\n")
cat(
  "Excluded chromosomes:",
  paste(excluded_chr, collapse = ", "),
  "\n"
)

if (length(prefixes) != 20) {
  warning(
    "Expected 20 autosomal chromosome files, but found ",
    length(prefixes)
  )
}

cat("\nFiles included in sample QC:\n")
print(basename(prefixes))

# ============================
# 3. Function to read one chromosome
# ============================

read_chr_stats <- function(prefix) {

  het_file <- paste0(prefix, ".het")
  imiss_file <- paste0(prefix, ".imiss")
  idepth_file <- paste0(prefix, ".idepth")

  required_files <- c(
    het_file,
    imiss_file,
    idepth_file
  )

  missing_files <- required_files[!file.exists(required_files)]

  if (length(missing_files) > 0) {
    stop(
      "Missing files for prefix:\n",
      prefix,
      "\n",
      paste(missing_files, collapse = "\n")
    )
  }

  # ----- Heterozygosity -----
  het <- read_table(
    het_file,
    show_col_types = FALSE
  ) %>%
    rename_with(~ gsub("\\(|\\)", "", .x)) %>%
    rename(
      HET_SITES = N_SITES
    ) %>%
    select(
      INDV,
      OHOM,
      HET_SITES
    )

  # ----- Missingness -----
  imiss <- read_table(
    imiss_file,
    show_col_types = FALSE
  ) %>%
    rename_with(~ gsub("\\(|\\)", "", .x)) %>%
    select(
      INDV,
      N_MISS,
      N_DATA
    )

  # ----- Mean depth -----
  idepth <- read_table(
    idepth_file,
    show_col_types = FALSE
  ) %>%
    rename_with(~ gsub("\\(|\\)", "", .x)) %>%
    rename(
      DEPTH_SITES = N_SITES
    ) %>%
    select(
      INDV,
      DEPTH_SITES,
      MEAN_DEPTH
    )

  # Merge the three tables
  reduce(
    list(het, imiss, idepth),
    full_join,
    by = "INDV"
  )
}

# ============================
# 4. Read all autosomes
# ============================

all_chr_data <- prefixes %>%
  map_df(read_chr_stats, .id = "chr_file")

cat(
  "\nNumber of chromosome-sample rows:",
  nrow(all_chr_data),
  "\n"
)

cat(
  "Number of unique samples:",
  n_distinct(all_chr_data$INDV),
  "\n"
)

# ============================
# 5. Collapse chromosomes
#    to one row per sample
# ============================

combined <- all_chr_data %>%
  group_by(INDV) %>%
  summarise(
    total_OHOM = sum(
      OHOM,
      na.rm = TRUE
    ),

    total_sites = sum(
      HET_SITES,
      na.rm = TRUE
    ),

    total_miss = sum(
      N_MISS,
      na.rm = TRUE
    ),

    total_data = sum(
      N_DATA,
      na.rm = TRUE
    ),

    depth_weighted_sum = sum(
      MEAN_DEPTH * DEPTH_SITES,
      na.rm = TRUE
    ),

    total_depth_sites = sum(
      DEPTH_SITES,
      na.rm = TRUE
    ),

    .groups = "drop"
  ) %>%
  mutate(
    missing_rate = if_else(
      total_data > 0,
      total_miss / total_data,
      NA_real_
    ),

    mean_depth = if_else(
      total_depth_sites > 0,
      depth_weighted_sum / total_depth_sites,
      NA_real_
    ),

    # Observed proportion of heterozygous sites
    prop_het = if_else(
      total_sites > 0,
      1 - total_OHOM / total_sites,
      NA_real_
    )
  ) %>%
  select(
    -depth_weighted_sum,
    -total_depth_sites
  )

cat(
  "Samples in combined autosomal QC table:",
  nrow(combined),
  "\n"
)

# ============================
# 6. Load metadata and merge
# ============================

meta <- read_csv(
  meta_file,
  show_col_types = FALSE
)

if (!sample_col %in% names(meta)) {
  stop(
    "Sample column '",
    sample_col,
    "' was not found in metadata."
  )
}

if (!group_var %in% names(meta)) {
  stop(
    "Grouping column '",
    group_var,
    "' was not found in metadata."
  )
}

if (anyDuplicated(meta[[sample_col]]) > 0) {
  stop("Duplicated sample IDs found in metadata.")
}

idf <- meta %>%
  inner_join(
    combined,
    by = setNames("INDV", sample_col)
  )

cat("Samples in metadata:", nrow(meta), "\n")
cat("Samples after metadata-QC merge:", nrow(idf), "\n")

if (nrow(idf) == 0) {
  stop("No sample IDs matched between metadata and QC results.")
}

# Check unmatched samples
missing_from_qc <- setdiff(
  meta[[sample_col]],
  combined$INDV
)

missing_from_meta <- setdiff(
  combined$INDV,
  meta[[sample_col]]
)

cat(
  "Metadata samples absent from QC:",
  length(missing_from_qc),
  "\n"
)

cat(
  "QC samples absent from metadata:",
  length(missing_from_meta),
  "\n"
)

# ============================
# 7. Calculate excess heterozygosity
# ============================

idf <- idf %>%
  mutate(
    qc_group = factor(.data[[group_var]])
  ) %>%
  group_by(qc_group) %>%
  mutate(
    group_mean_het = mean(
      prop_het,
      na.rm = TRUE
    ),

    excess_het = prop_het - group_mean_het
  ) %>%
  ungroup() %>%
  mutate(
    missing_fail = (
      is.na(missing_rate) |
        missing_rate >= missing_cutoff
    ),

    excess_het_fail = (
      is.na(excess_het) |
        excess_het > excess_het_cutoff
    ),

    HQ_sample = (
      !missing_fail &
        !excess_het_fail
    )
  )

# ============================
# 8. QC plots
# ============================

p_depth <- ggplot(
  idf,
  aes(
    x = qc_group,
    y = mean_depth
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    linewidth = 0.8
  ) +
  geom_jitter(
    aes(color = qc_group),
    width = 0.25,
    height = 0,
    size = 2,
    alpha = 0.8
  ) +
  scale_color_viridis_d(option = "D") +
  labs(
    title = "Mean Sequencing Depth per Sample",
    x = "Sex",
    y = "Mean sequencing depth"
  ) +
  theme_classic(base_size = 18) +
  theme(
    plot.title = element_text(
      size = 22,
      face = "bold",
      hjust = 0.5
    ),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.position = "none"
  )

p_missing <- ggplot(
  idf,
  aes(
    x = qc_group,
    y = missing_rate
  )
) +
  geom_hline(
    yintercept = missing_cutoff,
    linetype = "dashed",
    linewidth = 0.8
  ) +
  geom_boxplot(
    outlier.shape = NA,
    linewidth = 0.8
  ) +
  geom_jitter(
    aes(color = qc_group),
    width = 0.25,
    height = 0,
    size = 2,
    alpha = 0.8
  ) +
  scale_color_viridis_d(option = "D") +
  labs(
    title = "Missing Genotype Rate per Sample",
    x = "Sex",
    y = "Missing genotype rate"
  ) +
  theme_classic(base_size = 18) +
  theme(
    plot.title = element_text(
      size = 22,
      face = "bold",
      hjust = 0.5
    ),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.position = "none"
  )

p_het <- ggplot(
  idf,
  aes(
    x = qc_group,
    y = prop_het
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    linewidth = 0.8
  ) +
  geom_jitter(
    aes(color = qc_group),
    width = 0.25,
    height = 0,
    size = 2,
    alpha = 0.8
  ) +
  scale_color_viridis_d(option = "D") +
  labs(
    title = "Observed Autosomal Heterozygosity",
    x = "Sex",
    y = "Observed heterozygosity"
  ) +
  theme_classic(base_size = 18) +
  theme(
    plot.title = element_text(
      size = 22,
      face = "bold",
      hjust = 0.5
    ),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.position = "none"
  )

p_excess_het <- ggplot(
  idf,
  aes(
    x = qc_group,
    y = excess_het
  )
) +
  geom_hline(
    yintercept = excess_het_cutoff,
    linetype = "dashed",
    linewidth = 0.8
  ) +
  geom_boxplot(
    outlier.shape = NA,
    linewidth = 0.8
  ) +
  geom_jitter(
    aes(color = qc_group),
    width = 0.25,
    height = 0,
    size = 2,
    alpha = 0.8
  ) +
  scale_color_viridis_d(option = "D") +
  labs(
    title = "Excess Autosomal Heterozygosity",
    x = "Sex",
    y = "Heterozygosity minus sex-group mean"
  ) +
  theme_classic(base_size = 18) +
  theme(
    plot.title = element_text(
      size = 22,
      face = "bold",
      hjust = 0.5
    ),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    legend.position = "none"
  )

p_combined <- (
  p_depth + p_missing
) /
  (
    p_het + p_excess_het
  ) +
  plot_layout(guides = "collect")

# ============================
# 9. Save plots
# ============================

ggsave(
  filename = file.path(
    stat_dir,
    "QTL_sample_QC_sex_corrected_autosomes.png"
  ),
  plot = p_combined,
  width = 16,
  height = 12,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = file.path(
    stat_dir,
    "QTL_depth_missing_sex_corrected_autosomes.png"
  ),
  plot = p_depth + p_missing,
  width = 14,
  height = 7,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = file.path(
    stat_dir,
    "QTL_heterozygosity_sex_corrected_autosomes.png"
  ),
  plot = p_het + p_excess_het,
  width = 14,
  height = 7,
  dpi = 300,
  bg = "white"
)

# ============================
# 10. Save complete QC metrics
# ============================

write_csv(
  idf,
  file.path(
    stat_dir,
    "QTL_sample_QC_metrics_sex_corrected_autosomes.csv"
  )
)

# ============================
# 11. High-quality sample selection
# ============================

HQsamples <- idf %>%
  filter(HQ_sample) %>%
  pull(all_of(sample_col))

writeLines(
  HQsamples,
  con = file.path(
    stat_dir,
    "QTL_HQ_samples_sex_corrected_autosomes.txt"
  )
)

# Save excluded sample details
excluded_samples <- idf %>%
  filter(!HQ_sample) %>%
  select(
    all_of(sample_col),
    all_of(group_var),
    mean_depth,
    missing_rate,
    prop_het,
    excess_het,
    missing_fail,
    excess_het_fail
  )

write_csv(
  excluded_samples,
  file.path(
    stat_dir,
    "QTL_excluded_samples_sex_corrected_autosomes.csv"
  )
)

# ============================
# 12. Print QC summary
# ============================

cat("\n====================================\n")
cat("Sex-corrected autosomal QC completed\n")
cat("====================================\n")

cat(
  "Total samples analysed:",
  nrow(idf),
  "\n"
)

cat(
  "Samples failing missingness cutoff:",
  sum(idf$missing_fail, na.rm = TRUE),
  "\n"
)

cat(
  "Samples showing excess heterozygosity:",
  sum(idf$excess_het_fail, na.rm = TRUE),
  "\n"
)

cat(
  "Samples failing either criterion:",
  sum(!idf$HQ_sample, na.rm = TRUE),
  "\n"
)

cat(
  "High-quality samples retained:",
  length(HQsamples),
  "\n"
)

cat(
  "Missingness cutoff:",
  missing_cutoff,
  "\n"
)

cat(
  "Excess heterozygosity cutoff:",
  excess_het_cutoff,
  "\n"
)

cat(
  "Results saved in:",
  stat_dir,
  "\n"
)