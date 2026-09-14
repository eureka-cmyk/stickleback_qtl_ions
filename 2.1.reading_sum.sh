#!/bin/bash
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4g
#SBATCH --time=01:00:00
#SBATCH --job-name=multi_bam_qc
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/multi-bam-qc-%j.out

# 如果已加载了通用 Java 模块导致冲突，可先卸载
module unload Java 2>/dev/null || true
# 若系统里找不到 java，再加载集群提供的 Java
if ! command -v java >/dev/null 2>&1; then
  module load java-uoneasy/17.0.6
fi

# 你的 Qualimap 安装路径
QUALIMAP=~/app/qualimap_v2.3/qualimap

# 输入/输出目录
IN_DIR="/gpfs01/home/shxyt1/data/sticklebacks/bams/stickleback_v5.0.1_assembly/bamstats/QC/raw_bams"
OUT_DIR="${IN_DIR}/Multi-Bam-QC"
mkdir -p "$OUT_DIR"

SUMMARY="${OUT_DIR}/global_raw_report_custom.txt"
TMP_LIST="${OUT_DIR}/qualimap.tmp.txt"

echo "Input QC folder:  $IN_DIR"
echo "Output folder:    $OUT_DIR"

# 表头（与后续 R 脚本一致）
echo -e "sample\tbam_file\treads\tmapped_reads\tpercentage_mapped\tmn_coverage\tstd_coverage\tAve_map_qc\tdupl_reads" > "$SUMMARY"

# 遍历每个样本文件夹，读取 genome_results.txt 里关键信息
for d in "${IN_DIR}"/QTL*; do
  sample=$(basename "$d")
  qc="${d}/genome_results.txt"
  echo "Processing: ${sample}"
  if [[ ! -f "$qc" ]]; then
    echo "  (skip) QC file not found: $qc"
    continue
  fi

  bam_file=$(grep -m1 "bam file =" "$qc" | awk -F' = ' '{print $2}')
  reads=$(grep -m1 "number of reads =" "$qc" | awk -F' = ' '{print $2}')
  mapped=$(grep -m1 "number of mapped reads =" "$qc" | awk -F' = ' '{print $2}')
  mn_cov=$(grep -m1 "mean coverageData =" "$qc" | awk -F' = ' '{print $2}')
  std_cov=$(grep -m1 "std coverageData =" "$qc" | awk -F' = ' '{print $2}')
  map_q=$(grep -m1 "mean mapping quality =" "$qc" | awk -F' = ' '{print $2}')
  dup_reads=$(grep -m1 "number of duplicated reads (estimated) =" "$qc" | awk -F' = ' '{print $2}')

  # 百分比映射（去逗号后计算）
  reads_num=$(echo "$reads" | tr -d ',')
  mapped_num=$(echo "$mapped" | tr -d ',')
  if [[ -n "$reads_num" && "$reads_num" -gt 0 ]]; then
    pct=$(awk -v a="$mapped_num" -v b="$reads_num" 'BEGIN{printf "%.2f", (a/b*100)}')
  else
    pct="NA"
  fi

  echo -e "${sample}\t${bam_file}\t${reads}\t${mapped}\t${pct}\t${mn_cov}\t${std_cov}\t${map_q}\t${dup_reads}" >> "$SUMMARY"
done

echo "Summary table created: $SUMMARY"

# 生成 multi-bamqc 的输入表（样本→该样本 QC 目录）
awk -v OFS='\t' -v inpath="$IN_DIR" 'NR>1 {print $1, inpath"/"$1"/"}' "$SUMMARY" > "$TMP_LIST"

# 生成综合 HTML 报告（不指定 --feature-file，避免 bed 依赖）
"$QUALIMAP" multi-bamqc -d "$TMP_LIST" -outformat HTML -outdir "$OUT_DIR" -outfile global_raw_report.html

echo "Multi-Bam-QC HTML report: ${OUT_DIR}/global_raw_report.html"
