## Load packages
.libPaths(c("/gpfs01/home/shxyt1/R/x86_64-pc-linux-gnu-library/4.4", .libPaths()))
library(ggplot2)
library(patchwork)
library(ggrepel)

## Get command line argument
args <- commandArgs(trailingOnly = TRUE)

bam_QC <- read.table("/gpfs01/home/shxyt1/data/sticklebacks/bams/stickleback_v5.0.1_assembly/bamstats/QC/raw_bams/Multi-Bam-QC/global_raw_report_custom.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Set working directory
setwd("/gpfs01/home/shxyt1/data/sticklebacks/bams/stickleback_v5.0.1_assembly/bamstats/QC/raw_bams/Multi-Bam-QC")

## Read in sample metadata
sample_data <- read.csv("/gpfs01/home/shxyt1/data/sticklebacks/pheno/bigdata_Christophe_header_2025-04-28.csv")

## Merge data with sample data
bam_QC <- merge(bam_QC, sample_data[, c(1, 6:15)], by.x = "sample", by.y = "individual", all.x = TRUE)

## Convert character columns to numeric
bam_QC$mn_coverage   <- as.numeric(gsub("X", "", x = bam_QC$mn_coverage))
bam_QC$std_coverage  <- as.numeric(gsub("X", "", x = bam_QC$std_coverage))
bam_QC$reads         <- as.numeric(gsub(",", "", x = bam_QC$reads))
bam_QC$mapped_reads <- as.numeric(gsub(",", "", gsub("\\s*\\(.*\\)", "", bam_QC$mapped_reads)))
bam_QC$dupl_reads    <- as.numeric(gsub(",", "", x = bam_QC$dupl_reads))   

## Plot 1: overview
# p <- ggplot(bam_QC) +
#   geom_point(aes(mn_coverage, mapped_reads, col = Ecotype))
# q <- ggplot(bam_QC) +
#   geom_boxplot(aes(Ecotype, mn_coverage, col = Ecotype), outlier.shape = NA) +
#   geom_jitter(aes(Ecotype, mn_coverage, col = Ecotype), height = 0) +
#   geom_hline(yintercept = 5, col = "red")
# s <- ggplot(bam_QC) +
#   geom_boxplot(aes(Ecotype, Ave_map_qc, col = Ecotype), outlier.shape = NA) +
#   geom_jitter(aes(Ecotype, Ave_map_qc, col = Ecotype), height = 0) +
#   geom_hline(yintercept = 10, col = "red")
# r <- ggplot(bam_QC) +
#   geom_boxplot(aes(Population, mn_coverage, col = Ecotype), outlier.shape = NA) +
#   geom_jitter(aes(Population, mn_coverage, col = Ecotype), height = 0) +
#   geom_hline(yintercept = 5, col = "red") +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
# t <- ggplot(bam_QC) +
#   geom_boxplot(aes(Ecotype, dupl_reads / mapped_reads, col = Ecotype), outlier.shape = NA) +
#   geom_jitter(aes(Ecotype, dupl_reads / mapped_reads, col = Ecotype), height = 0)


min_cov <- 5
min_QC <- 10
min_mapped_reads <- 10000

bam_QC$failed_sample <- with(
  bam_QC,
  mn_coverage < min_cov | Ave_map_qc < min_QC | mapped_reads < min_mapped_reads
)


p <- ggplot(bam_QC) +
  geom_point(aes(mn_coverage, mapped_reads, col = Ecotype, shape = failed_sample)) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4)) +  # TRUE=×, FALSE=圆点
  labs(title = "Mean Coverage vs Mapped Reads") +
  theme_bw()

q <- ggplot(bam_QC) +
  geom_boxplot(aes(Ecotype, mn_coverage, col = Ecotype), outlier.shape = NA) +
  geom_jitter(aes(Ecotype, mn_coverage, col = Ecotype, shape = failed_sample), height = 0) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4)) +
  geom_hline(yintercept = 5, col = "red") +
  labs(title = "Mean Coverage by Ecotype") +
  theme_bw()

s <- ggplot(bam_QC) +
  geom_boxplot(aes(Ecotype, Ave_map_qc, col = Ecotype), outlier.shape = NA) +
  geom_jitter(aes(Ecotype, Ave_map_qc, col = Ecotype, shape = failed_sample), height = 0) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4)) +
  geom_hline(yintercept = 10, col = "red") +
  labs(title = "Mean Mapping Quality by Ecotype") +
  theme_bw()

r <- ggplot(bam_QC) +
  geom_boxplot(aes(Population, mn_coverage, col = Ecotype), outlier.shape = NA) +
  geom_jitter(aes(Population, mn_coverage, col = Ecotype, shape = failed_sample), height = 0) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4)) +
  geom_hline(yintercept = 5, col = "red") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(title = "Mean Coverage by Population") +
  theme_bw()

t <- ggplot(bam_QC) +
  geom_boxplot(aes(Ecotype, dupl_reads / mapped_reads, col = Ecotype), outlier.shape = NA) +
  geom_jitter(aes(Ecotype, dupl_reads / mapped_reads, col = Ecotype, shape = failed_sample), height = 0) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4)) +
  labs(title = "Duplicate Reads Ratio by Ecotype") +
  theme_bw()

p1 <- (p / q / s / r / t)
ggsave(filename = "Mapping_coverage_mark.png", p1, width = 10, height = 20)


### Filter out low-quality samples
min_cov <- 5
min_QC <- 10
min_mapped_reads <- 10000

print(paste0(sum(as.numeric(!bam_QC$mn_coverage >= min_cov)), " samples have mean coverage < ", min_cov, "X"))
print(paste0(sum(as.numeric(!bam_QC$Ave_map_qc >= min_QC)), " samples have mean mapping QC < ", min_QC))
print(paste0("Removing ", sum(as.numeric(!bam_QC$mn_coverage >= min_cov & bam_QC$Ave_map_qc >= min_QC & bam_QC$mapped_reads >= min_mapped_reads)),
             " samples; retaining ", sum(as.numeric(bam_QC$mn_coverage >= min_cov & bam_QC$Ave_map_qc >= min_QC & bam_QC$mapped_reads >= min_mapped_reads))))

bam_QC_hiQ <- bam_QC[bam_QC$mn_coverage >= min_cov & bam_QC$Ave_map_qc >= min_QC & bam_QC$mapped_reads >= min_mapped_reads, ]

print(paste0("Mean mapped reads = ", mean(bam_QC_hiQ$mapped_reads),
             " (min=", min(bam_QC_hiQ$mapped_reads), ", max=", max(bam_QC_hiQ$mapped_reads), ")"))
print(paste0("Mean coverage = ", mean(bam_QC_hiQ$mn_coverage),
             " (min=", min(bam_QC_hiQ$mn_coverage), ", max=", max(bam_QC_hiQ$mn_coverage), ")"))
print(paste0("Mean mapping QC = ", mean(bam_QC_hiQ$Ave_map_qc),
             " (min=", min(bam_QC_hiQ$Ave_map_qc), ", max=", max(bam_QC_hiQ$Ave_map_qc), ")"))



#bam_QC_hiQ
p <- ggplot(bam_QC_hiQ) +
  geom_point(aes(mn_coverage, mapped_reads, col = Ecotype))
q <- ggplot(bam_QC_hiQ) +
  geom_boxplot(aes(Ecotype, mn_coverage, col = Ecotype), outlier.shape = NA) +
  geom_jitter(aes(Ecotype, mn_coverage, col = Ecotype), height = 0) +
  geom_hline(yintercept = 5, col = "red")
s <- ggplot(bam_QC_hiQ) +
  geom_boxplot(aes(Ecotype, Ave_map_qc, col = Ecotype), outlier.shape = NA) +
  geom_jitter(aes(Ecotype, Ave_map_qc, col = Ecotype), height = 0) +
  geom_hline(yintercept = 10, col = "red")
r <- ggplot(bam_QC_hiQ) +
  geom_boxplot(aes(Population, mn_coverage, col = Ecotype), outlier.shape = NA) +
  geom_jitter(aes(Population, mn_coverage, col = Ecotype), height = 0) +
  geom_hline(yintercept = 5, col = "red") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
t <- ggplot(bam_QC_hiQ) +
  geom_boxplot(aes(Ecotype, dupl_reads / mapped_reads, col = Ecotype), outlier.shape = NA) +
  geom_jitter(aes(Ecotype, dupl_reads / mapped_reads, col = Ecotype), height = 0)

p1 <- (p / q / s / r / t)
ggsave(filename = "Mapping_coverage_HiQ.png", p1, width = 10, height = 20)




# Write output files
write.table(bam_QC_hiQ, file = "HiQ_bam_files.csv",
            row.names = FALSE, quote = FALSE, col.names = TRUE, sep = ",")
writeLines(con = "HiQ_bam_files.txt", paste(bam_QC_hiQ$sample, bam_QC_hiQ$bam_file))
