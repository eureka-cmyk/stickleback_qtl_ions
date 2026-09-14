#!/bin/bash
#SBATCH --job-name=vcf_qc
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/qc-%A_%a.out
#SBATCH --array=1-22

module load vcftools-uoneasy/0.1.16-GCC-12.3.0

# === 输入与路径设置 ===
vcf_dir="/share/MacColl_shared/QTL/vcf/GCF_016920845.1_GAculeatus_UGA_version5_genomic/ploidy_aware_HWEPops_MQ10_BQ20_with_sex"
stat_dir="/share/MacColl_shared/QTL/vcf/stats_with_sex"

mkdir -p $stat_dir

# === 获取当前 array 的第 N 条染色体 BCF ===
# 自动按顺序列出所有 *_sorted.bcf 并取第 $SLURM_ARRAY_TASK_ID 个
bcf_file=$(ls $vcf_dir/*_sorted.bcf | sed -n "${SLURM_ARRAY_TASK_ID}p")

if [[ ! -f "$bcf_file" ]]; then
    echo "Error: No BCF found for array index $SLURM_ARRAY_TASK_ID"
    exit 1
fi

echo "Processing: $bcf_file"

# 获取文件基名，不含 .bcf
base=$(basename $bcf_file .bcf)
stat_prefix="${stat_dir}/${base}"

# === 运行 VCFtools QC ===

echo "Running depth..."
vcftools --bcf $bcf_file --depth --out ${stat_prefix} &

echo "Running missing-indv..."
vcftools --bcf $bcf_file --missing-indv --out ${stat_prefix} &

echo "Running het..."
vcftools --bcf $bcf_file --het --out ${stat_prefix} &

echo "Running relatedness..."
vcftools --bcf $bcf_file --relatedness --out ${stat_prefix} &

# 等待所有后台进程完成
wait

echo "Finished QC for $bcf_file"
