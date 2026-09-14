library(data.table)
library(CMplot)

files <- list(
  Ca = "/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_nine_main/Ca_no_outlier/Ca_no_outlier_MLMA.mlma",
  Co = "/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_nine_main/Co_no_outlier/Co_no_outlier_MLMA.mlma",
  Fe = "/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_nine_main/Fe_no_outlier/Fe_no_outlier_MLMA.mlma",
  Sr = "/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_nine_main/Sr_no_outlier/Sr_no_outlier_MLMA.mlma",
  P  = "/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_no_outlier/P_log/P_log_MLMA.mlma",
  V  = "/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_no_outlier/V_log/V_log_MLMA.mlma"
)

read_mlma <- function(f) {
  dt <- fread(f)
  # 标准化到 CMplot 需要的列名
  dt2 <- dt[, .(
    SNP = as.character(SNP),
    CHR = as.integer(Chr),
    POS = as.integer(bp),
    P   = as.numeric(p)
  )]
  dt2 <- dt2[!is.na(CHR) & !is.na(POS) & !is.na(P)]
  dt2
}

# 读入并把每个性状的 P 列改成 trait 名
dt_list <- lapply(names(files), function(tr) {
  dt <- read_mlma(files[[tr]])
  setnames(dt, "P", tr)   # P -> Ca/Co/...
  dt
})
names(dt_list) <- names(files)

# 用 data.table 的 merge 合并（全外连接）
merged <- Reduce(function(x, y) merge(x, y, by = c("SNP","CHR","POS"), all = TRUE),
                 dt_list)

setorder(merged, CHR, POS)

# 画环状曼哈顿（6圈）
# CMplot(
#   merged,
#   plot.type="c",
#   multracks=TRUE,
#   LOG10=TRUE,
#   threshold=c(5e-8),
#   threshold.col="red",
#   r=0.45,
#   cir.chr.h=1.3,
#   cir.axis=TRUE,
#   outward=FALSE,
#   chr.den.col="black",
#   file="pdf",
#   file.name="Stickleback_6traits_circos_manhattan",
#   dpi=300,
#   width=10,
#   height=10,
#   verbose=TRUE
# )

CMplot(
  merged,
  plot.type   = "c",      # 环状 Manhattan
  multracks   = TRUE,     # 6 个性状 = 6 个环
  LOG10       = TRUE,
  
  threshold       = c(5e-8),
  threshold.col   = "red",
  
  # ---------- 点的控制（关键：瘦身） ----------
  cex            = 0.25,        # 普通 SNP 点更小
  signal.cex     = 0.5,         # 显著 SNP 稍大
  signal.col     = rep("red", 6),  # 所有性状显著点标红
  
  # ---------- 圆环参数 ----------
  r              = 0.42,        # 中心空白略大 → 减少拥挤
  cir.chr.h      = 1.15,        # 环间距略小
  cir.axis       = TRUE,
  outward        = FALSE,
  chr.den.col    = "black",
  
  # ---------- 输出 ----------
  file           = "png",       # <<<<<< 关键：不要再用 pdf
  file.name      = "Stickleback_6traits_circos_manhattan",
  dpi            = 200,         # 足够清晰，文件体积明显下降
  width          = 8,           # 英寸
  height         = 8,
  verbose        = TRUE
)

