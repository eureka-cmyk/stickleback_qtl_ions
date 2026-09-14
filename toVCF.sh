#!/bin/bash
# Yiwen
# 4/11/25
# Variant calling for all BAMs on Ada HPC

#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64g
#SBATCH --time=72:00:00
#SBATCH --job-name=bwa_mapping
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/toVCF-%x-%j.out

# 加载模块
module load BCFtools/1.18-GCC-12.3.0
module load samtools-uoneasy/1.18-GCC-12.3.0

# 路径设置
bam_dir="/gpfs01/home/shxyt1/data/sticklebacks/bams/stickleback_v5.0.1_assembly/raw_bams"
ref_genome="/gpfs01/home/shxyt1/data/sticklebacks/ref/fresh/stickleback_v5.0.1_assembly.fa"
out_dir="/gpfs01/home/shxyt1/data/sticklebacks/VCF/stickleback_v5.0.1_assembly"

mkdir -p "$out_dir"

# 输出文件名
vcf_prefix="${out_dir}/stickleback_all"

echo "Running variant calling for all BAMs in ${bam_dir}"
echo "Reference genome: ${ref_genome}"
echo "Output prefix: ${vcf_prefix}"

# Step 1: 生成变异（mpileup + call）
bcftools mpileup -Ou \
    -f "$ref_genome" \
    -b <(ls "$bam_dir"/*.bam) \
    --threads $SLURM_CPUS_PER_TASK \
    | bcftools call -mv -Oz -o "${vcf_prefix}.raw.vcf.gz"

# Step 2: 索引
bcftools index "${vcf_prefix}.raw.vcf.gz"

# Step 3: 基础过滤（可选，根据需求修改）
bcftools view -i 'QUAL>30 && DP>10' "${vcf_prefix}.raw.vcf.gz" -Oz -o "${vcf_prefix}.filtered.vcf.gz"
bcftools index "${vcf_prefix}.filtered.vcf.gz"

echo "✅ Variant calling complete!"
echo "Output files:"
echo "  - ${vcf_prefix}.raw.vcf.gz"
echo "  - ${vcf_prefix}.filtered.vcf.gz"
