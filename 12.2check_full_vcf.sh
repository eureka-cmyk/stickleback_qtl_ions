#!/bin/bash
#SBATCH --job-name=check_merge
#SBATCH --partition=defq
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/check_merge_%j.log

#set -euo pipefail

# =========================
# Load required module
# =========================
module load bcftools-uoneasy/1.19-GCC-13.2.0

# =========================
# Explicit path to final merged VCF
# =========================
VCF_FULL="/gpfs01/home/shxyt1/data/sticklebacks/VCF/merge/stickleback_allchr_full.vcf.gz"

echo "===== [1] Check file exists ====="
ls -lh "${VCF_FULL}"

echo "===== [2] Force re-index VCF (tabix) ====="
tabix -f -p vcf "${VCF_FULL}"

echo "===== [3] Sample count (should be 352) ====="
bcftools query -l "${VCF_FULL}" | wc -l

echo "===== [4] Chromosome list (FULL, no truncation) ====="
bcftools view -H "${VCF_FULL}" \
  | cut -f1 \
  | sort -u

echo "===== [5] bcftools stats (header preview) ====="
bcftools stats "${VCF_FULL}" | head -n 40

echo "===== CHECK COMPLETE ====="
