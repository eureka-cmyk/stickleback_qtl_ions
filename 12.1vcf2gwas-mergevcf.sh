#!/bin/bash
#SBATCH --job-name=merge_manhattan
#SBATCH --partition=defq
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/merge_manhattan_%j.log


module load bcftools-uoneasy/1.19-GCC-13.2.0

#bcftools concat \
#/gpfs01/home/shxyt1/data/sticklebacks/VCF/stickleback_NC_05321*.vcf.gz \
#-Oz \
#-o /gpfs01/home/shxyt1/data/sticklebacks/VCF/merge/stickleback_allchr.vcf.gz

VCF_DIR=/gpfs01/home/shxyt1/data/sticklebacks/VCF
MERGE_DIR=/gpfs01/home/shxyt1/data/sticklebacks/VCF/merge


bcftools concat \
${MERGE_DIR}/stickleback_allchr.vcf.gz \
${VCF_DIR}/stickleback_NC_05322*.vcf.gz \
${VCF_DIR}/stickleback_NC_05323*.vcf.gz \
-Oz \
-o ${MERGE_DIR}/stickleback_allchr_full.vcf.gz




tabix -p vcf \
/gpfs01/home/shxyt1/data/sticklebacks/VCF/merge/stickleback_allchr_full.vcf.gz


# 样本数
bcftools query -l /gpfs01/home/shxyt1/data/sticklebacks/VCF/merge/stickleback_allchr_full.vcf.gz | wc -l

# 染色体是否齐全
bcftools view -H /gpfs01/home/shxyt1/data/sticklebacks/VCF/merge/stickleback_allchr_full.vcf.gz| cut -f1 | sort | uniq | head

# 基本统计
bcftools stats /gpfs01/home/shxyt1/data/sticklebacks/VCF/merge/stickleback_allchr_full.vcf.gz | head



