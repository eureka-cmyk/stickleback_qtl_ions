.libPaths(c("/gpfs01/home/shxyt1/R/x86_64-pc-linux-gnu-library/4.4", .libPaths()))
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(patchwork)

inputdir <- "/share/MacColl_shared/QTL/bams/GCF_016920845.1_GAculeatus_UGA_version5_genomic/Sex_coverage/"

# coverage file
files <- list.files(inputdir, pattern = "_coverage.txt$", full.names = TRUE)

# read and combine
cov <- files %>%
  map_df(~ read.table(.x, header = FALSE, sep = "\t",
                      col.names = c("individual", "chrom", "value")))

head(cov)
# trans
cov_wide <- cov %>%
  pivot_wider(names_from = chrom, values_from = value)

# ratio
cov_wide <- cov_wide %>%
  mutate(
    X_ratio = NC_053230.1 / NC_053212.1,
    Y_ratio = NC_053233.1 / NC_053212.1
  )

# threshold
X_threshold <- 0.85

cov_wide <- cov_wide %>%
  mutate(
    sex = ifelse(X_ratio >= X_threshold, "F", "M")
  )

# save
write.table(cov_wide,
            file = file.path(inputdir, "Sex_calling_summary.csv"),
            sep = ",", quote = FALSE, row.names = FALSE)

# PED
ped <- cov_wide %>%
  mutate(fam="Fam1", dad=0, mom=0) %>%
  select(fam, individual, dad, mom, sex)

write.table(ped,
            file = file.path(inputdir, "Gsex.ped"),
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

# ===== pkot =====

p1 <- ggplot(cov_wide, aes(individual, X_ratio, color = sex)) +
  geom_point(size = 3) +
  geom_hline(yintercept = X_threshold, linetype = "dashed", color="red") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust=1)) +
  ggtitle("X Coverage Ratio (X/Autosome)")

p2 <- ggplot(cov_wide, aes(individual, Y_ratio, color = sex)) +
  geom_point(size = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust=1)) +
  ggtitle("Y Coverage Ratio (Y/Autosome)")

p3 <- ggplot(cov_wide, aes(X_ratio, Y_ratio, color = sex, label = individual)) +
  geom_point(size = 3) +
  geom_text_repel() +
  theme_bw() +
  ggtitle("Sex Identification: X_ratio vs Y_ratio")

# save plot
ggsave(file.path(inputdir, "Sex_cov_X_ratio.png"), p1, width=10, height=5)
ggsave(file.path(inputdir, "Sex_cov_Y_ratio.png"), p2, width=10, height=5)
ggsave(file.path(inputdir, "Sex_cov_XY_scatter.png"), p3, width=7, height=6)

# 
big_plot <- (p1 + p2) / p3
ggsave(file.path(inputdir, "Sex_cov_summary.png"), big_plot, width=12, height=10)


#QTL22197 QTL22254

ydir <- "/gpfs01/home/shxyt1/data/sticklebacks/pheno/Y_depth"

bin_size <- 50000

files <- list.files(ydir, pattern = "_NC053233_depth.txt$", full.names = TRUE)

ydepth <- map_df(files, function(f){
  id <- sub("_NC053233_depth.txt", "", basename(f))
  read.table(f, header = FALSE, col.names = c("chr","pos","depth")) %>%
    mutate(ID = id)
})

ybin <- ydepth %>%
  mutate(bin = (pos - 1) %/% bin_size) %>%
  group_by(ID, bin) %>%
  summarise(
    mean_depth = mean(depth, na.rm = TRUE),
    start_pos  = bin * bin_size + 1,
    .groups = "drop"
  ) %>%
  mutate(pos_Mb = start_pos / 1e6)

p <- ggplot(ybin, aes(x = pos_Mb, y = mean_depth, color = ID)) +
  geom_line() +
  labs(
    x = "Position on Y chromosome (Mb)",
    y = paste0("Mean depth per ", bin_size/1000, " kb window"),
    title = "Coverage across Y chromosome (NC_053233.1)"
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(file.path(ydir, "Y_coverage_QTL22197_QTL22254.png"),
       p, width = 9, height = 5)

