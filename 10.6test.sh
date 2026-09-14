#!/bin/bash
#SBATCH --job-name=gemma_test_5pct
#SBATCH --partition=defq
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/gemma_test_5pct_%j.log

echo "Start GEMMA 5% thin test"
date
hostname

module load plink-uoneasy/1.9b_6.21-x86_64

# ---------- PATHS ----------
PLINK="/gpfs01/home/shxyt1/data/sticklebacks/PLINK/stickleback_all_clean"
OUTDIR="/gpfs01/home/shxyt1/data/sticklebacks/gemma/Na5%"
PHENO_FILE="/gpfs01/home/shxyt1/data/sticklebacks/pheno/traits/pheno_Na_RINT.txt"

mkdir -p $OUTDIR

cd $OUTDIR
echo "Working in: $(pwd)"

# ------------------------------
# 1. THIN 5% SNPs (random subset)
# ------------------------------
echo "[1] Randomly thinning SNPs to 5%..."
plink --bfile $PLINK --thin 0.05 --make-bed --out stickleback_thin5

echo "[1] Done thinning."
ls -lh stickleback_thin5.*

# ------------------------------
# 2. Run GEMMA BSLMM
# ------------------------------
echo "[2] Running GEMMA BSLMM..."

gemma \
  -bfile stickleback_thin5 \
  -bslmm 1 \
  -o Na5pct_BSLMM \
  -p $PHENO_FILE

echo "[2] GEMMA finished."
date

echo "Results saved in:"
echo "  $OUTDIR/output/"
