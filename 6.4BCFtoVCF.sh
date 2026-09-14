#!/bin/bash
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --job-name=sexcov
#SBATCH --time=04:00:00
#SBATCH --mem=4G
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/BCF-vcf-%j.out

module load bcftools-uoneasy/1.19-GCC-13.2.0

bcf_dir=/share/MacColl_shared/QTL/vcf/GCF_016920845.1_GAculeatus_UGA_version5_genomic/ploidy_aware_HWEPops_MQ10_BQ20
out_dir=/share/MacColl_shared/QTL/vcf/GCF_016920845.1_GAculeatus_UGA_version5_genomic/sex_coverage
mkdir -p "$out_dir"

# 1. autosome NC_053212.1
bcftools view -O z -o "$out_dir/stickleback_NC_053212.1_sorted.vcf.gz"  \
  "$bcf_dir/stickleback_NC_053212.1_sorted.bcf"
bcftools index -t "$out_dir/stickleback_NC_053212.1_sorted.vcf.gz"

# 2. X NC_053230.1
bcftools view -O z -o "$out_dir/stickleback_NC_053230.1_sorted.vcf.gz"  \
  "$bcf_dir/stickleback_NC_053230.1_sorted.bcf"
bcftools index -t "$out_dir/stickleback_NC_053230.1_sorted.vcf.gz"

# 3. Y NC_053233.1
bcftools view -O z -o "$out_dir/stickleback_NC_053233.1_sorted.vcf.gz"  \
  "$bcf_dir/stickleback_NC_053233.1_sorted.bcf"
bcftools index -t "$out_dir/stickleback_NC_053233.1_sorted.vcf.gz"
