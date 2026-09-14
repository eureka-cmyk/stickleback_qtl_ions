library(tidyverse)
library(patchwork)
library(viridis)
install.packages("viridis")
# ============================
# 1. Define paths
# ============================

stat_dir <- "/share/MacColl_shared/QTL/vcf/stats_with_sex"
meta_file <- "/gpfs01/home/shxyt1/data/sticklebacks/pheno/qtl_all_phenotype_data.csv"
sample_col <- "fish_id"

# ============================
# 2. List chromosome prefixes
# ============================

het_files <- list.files(stat_dir, pattern="\\.het$", full.names=TRUE)
prefixes <- gsub("\\.het$", "", het_files)

# ============================
# 3. Function to read 1 chromosome stats
# ============================

read_chr_stats <- function(prefix) {
  
  # ----- Read het -----
  het <- read_table(paste0(prefix, ".het")) %>%
    rename_with(~ gsub("\\(|\\)", "", .x)) %>%
    rename(HET_SITES = N_SITES) %>%     # rename N_SITES
    select(INDV, OHOM, HET_SITES)
  
  # ----- Read missing -----
  imiss <- read_table(paste0(prefix, ".imiss")) %>%
    rename_with(~ gsub("\\(|\\)", "", .x)) %>%
    select(INDV, N_MISS, N_DATA)
  
  # ----- Read depth -----
  idepth <- read_table(paste0(prefix, ".idepth")) %>%
    rename_with(~ gsub("\\(|\\)", "", .x)) %>%
    rename(DEPTH_SITES = N_SITES) %>%   # rename N_SITES
    select(INDV, DEPTH_SITES, MEAN_DEPTH)
  
  # Merge 3 tables for this chromosome
  df <- reduce(list(het, imiss, idepth), merge, by="INDV")
  return(df)
}

# ============================
# 4. Read all chromosomes
# ============================

all_chr_data <- prefixes %>%
  map_df(read_chr_stats)

# ============================
# 5. Collapse 22 chromosomes -> 1 row per sample
# ============================

combined <- all_chr_data %>%
  group_by(INDV) %>%
  summarise(
    total_OHOM  = sum(OHOM),
    total_sites = sum(HET_SITES),
    total_miss  = sum(N_MISS),
    total_data  = sum(N_DATA),
    # Weighted depth correct calculation
    mean_depth  = sum(MEAN_DEPTH * DEPTH_SITES) / sum(DEPTH_SITES)
  ) %>%
  mutate(
    missing_rate = total_miss / total_data,
    prop_het = 1 - (total_OHOM / total_sites)
  )

# ============================
# 6. Load metadata & merge
# ============================

meta <- read.csv(meta_file, header=TRUE)
idf <- merge(meta, combined, by.x=sample_col, by.y="INDV")

# ============================
# 7. Group by sex (or population if available)
# ============================

group_var <- "sex"

idf <- idf %>%
  group_by(.data[[group_var]]) %>%
  mutate(
    excess_het = prop_het - mean(prop_het, na.rm=TRUE),
    excess_het_flag = excess_het > 0.015
  )

# ============================
# 8. QC plots
# ============================

p_depth <- ggplot(idf, aes(x = .data[[group_var]], y = mean_depth)) +
  geom_boxplot(outlier.shape = NA, size=1) +
  geom_jitter(aes(color = .data[[group_var]]),
              width=0.25, height=0, size=2, alpha=0.8) +
  scale_color_viridis_d(option="D") +
  theme_classic(base_size = 18) +
  ggtitle("Mean Sequencing Depth per Sample") +
  theme(
    plot.title = element_text(size = 22, face = "bold"),
    axis.text = element_text(size=16),
    axis.title = element_text(size=18),
    legend.position = "none"   # 可根据需要打开
  )


p_missing <- ggplot(idf, aes(x = .data[[group_var]], y = missing_rate)) +
  geom_boxplot(outlier.shape = NA, size=1) +
  geom_jitter(aes(color = .data[[group_var]]),
              width=0.25, height=0, size=2, alpha=0.8) +
  scale_color_viridis_d(option="D") +
  theme_classic(base_size = 18) +
  ggtitle("Missing Genotype Rate per Sample") +
  theme(
    plot.title = element_text(size = 22, face = "bold"),
    axis.text = element_text(size=16),
    axis.title = element_text(size=18),
    legend.position = "none"
  )



p_het <- ggplot(idf, aes(x = .data[[group_var]], y = prop_het)) +
  geom_boxplot(outlier.shape = NA, size=1) +
  geom_jitter(color = viridis(1, option="D"),
              width=0.25, height=0, size=2, alpha=0.8) +
  theme_classic(base_size = 18) +
  ggtitle("Observed Heterozygosity per Sample") +
  theme(
    plot.title = element_text(size = 22, face = "bold"),
    axis.text = element_text(size=16),
    axis.title = element_text(size=18)
  )



p_combined <- p_depth + p_missing + plot_layout(guides="collect")

# ============================
# 9. Save plots
# ============================

ggsave(paste0(stat_dir, "/QTL_depth_missing_corrected.png"),
       p_combined, width=14, height=10)

ggsave(paste0(stat_dir, "/QTL_heterozygosity_corrected.png"),
       p_het, width=14, height=10)

# ============================
# 10. High-quality sample selection
# ============================

HQsamples <- idf[[sample_col]][
  idf$missing_rate < 0.4 &
    idf$excess_het < 0.015
]

writeLines(HQsamples,
           con=paste0(stat_dir, "/QTL_HQ_samples_corrected.txt"))

cat(sum(idf$missing_rate >= 0.4), "samples fail missingness cutoff\n")
cat(sum(idf$excess_het > 0.015), "samples show excess heterozygosity\n")
#90 samples show excess heterozygosity
cat("Correct QC completed.\n")

