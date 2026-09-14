#!/bin/bash
#SBATCH --job-name=plink_merge_22chr
#SBATCH --partition=defq
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/plink_merge_22chr_%A.log

module load plink-uoneasy/1.9b_6.21-x86_64

PLINK_DIR="/gpfs01/home/shxyt1/data/sticklebacks/PLINK"
cd "${PLINK_DIR}"

echo "Working directory: $(pwd)"
echo "Listing per-chromosome PLINK files..."

# 找到所有 per-chr 前缀
ls stickleback_NC_*.bed | sed 's/\.bed$//' | sort > all_plink_prefixes.txt

echo "All chromosome prefixes:"
cat all_plink_prefixes.txt

# 第一条作为 base
BASE_PREFIX=$(head -n 1 all_plink_prefixes.txt)
tail -n +2 all_plink_prefixes.txt > merge_list.txt

echo "Base dataset    : ${BASE_PREFIX}"
echo "Merge-list file : merge_list.txt"
echo "Content of merge_list.txt:"
cat merge_list.txt

# 合并
plink \
  --bfile "${BASE_PREFIX}" \
  --merge-list merge_list.txt \
  --allow-extra-chr \
  --make-bed \
  --out stickleback_allchr

echo "Done. Merged dataset generated:"
echo "  ${PLINK_DIR}/stickleback_allchr.bed"
echo "  ${PLINK_DIR}/stickleback_allchr.bim"
echo "  ${PLINK_DIR}/stickleback_allchr.fam"
