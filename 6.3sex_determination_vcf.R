library(data.table)
library(dplyr)

dp_dir <- "/share/MacColl_shared/QTL/vcf/GCF_016920845.1_GAculeatus_UGA_version5_genomic/sex_coverage_DP"
mean_dir <- file.path(dp_dir, "mean")

# ----------------------------------------
# 1. 获取 VCF 原始 sample 顺序
# ----------------------------------------
vcf_file <- "/share/MacColl_shared/QTL/vcf/GCF_016920845.1_GAculeatus_UGA_version5_genomic/ploidy_aware_HWEPops_MQ10_BQ20/stickleback_NC_053212.1_sorted.bcf"
sample_list <- system(paste("bcftools query -l", vcf_file), intern = TRUE)

length(sample_list)  # 351

# ----------------------------------------
# 2. 计算 autosome / X / Y 的 mean depth
# ----------------------------------------
chr_auto <- "NC_053212.1"
chr_X <- "NC_053230.1"
chr_Y <- "NC_053233.1"

get_mean_depth <- function(chr){
  file_dp <- file.path(dp_dir, paste0(chr, "_DP.txt"))
  message("Reading ", file_dp)
  
  dt <- fread(file_dp)
  # 第 3 列以后是 DP, 对应 sample_list
  depth_mat <- as.matrix(dt[, 3:ncol(dt)])
  mean_depth <- colMeans(depth_mat, na.rm = TRUE)
  
  tibble(
    ID = sample_list,
    mean = mean_depth
  )
}

auto <- get_mean_depth(chr_auto)
Xcov <- get_mean_depth(chr_X)
Ycov <- get_mean_depth(chr_Y)

# ----------------------------------------
# 3. 合并成 ratio 表格
# ----------------------------------------
cov_df <- auto %>%
  rename(auto = mean) %>%
  left_join(Xcov %>% rename(X = mean), by="ID") %>%
  left_join(Ycov %>% rename(Y = mean), by="ID") %>%
  mutate(
    X_ratio = X / auto,
    Y_ratio = Y / auto
  )

write.csv(cov_df, file.path(dp_dir, "sex_cov_ratio_corrected.csv"), row.names = FALSE)

message("Finished! ✔ Sample 顺序已完全修复。")
