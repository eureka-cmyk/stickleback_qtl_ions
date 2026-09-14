#!/bin/bash
#SBATCH --job-name=plink_check
#SBATCH --partition=defq
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/plink_check_%j.log

module load plink-uoneasy/1.9b_6.21-x86_64

# ------------------------------------------------------------
#  INPUT AND OUTPUT PATHS
# ------------------------------------------------------------
PLINK_MASTER="/gpfs01/home/shxyt1/data/sticklebacks/PLINK/stickleback_allchr"
OUTDIR="/gpfs01/home/shxyt1/data/sticklebacks/plink_check"
mkdir -p $OUTDIR

echo "Starting GWAS preparation pipeline..."
echo "Input dataset: $PLINK_MASTER"
echo "Output folder: $OUTDIR"


# ------------------------------------------------------------
# 1. BASIC QC: FILTER INDIVIDUALS & SNPs BY MISSINGNESS
#    - Remove SNPs with missing rate > 5%
#    - Remove individuals with missing rate > 10%
# ------------------------------------------------------------
echo "[1] Running basic QC (missingness filters)..."

plink \
  --bfile $PLINK_MASTER \
  --geno 0.05 \
  --mind 0.10 \
  --make-bed \
  --out $OUTDIR/step1_basicQC


# ------------------------------------------------------------
# 2. SEX CHECK (Optional: only works if .fam contains sexes)
#    - Creates sexcheck file
#    - You can manually remove mismatches afterwards
# ------------------------------------------------------------
echo "[2] Performing sex check..."

plink \
  --bfile $OUTDIR/step1_basicQC \
  --check-sex \
  --out $OUTDIR/step2_sexcheck


# ------------------------------------------------------------
# 3. LD PRUNING
#    - Creates a set of independent SNPs for PCA
#    - Window size = 50 SNPs
#    - Step size = 5 SNPs
#    - r^2 threshold = 0.2
# ------------------------------------------------------------
echo "[3] Running LD pruning..."

plink \
  --bfile $OUTDIR/step1_basicQC \
  --indep-pairwise 50 5 0.2 \
  --out $OUTDIR/step3_prunedSNPs


# ------------------------------------------------------------
# 4. PCA USING PRUNED SNP SET
# ------------------------------------------------------------
echo "[4] Performing PCA..."

plink \
  --bfile $OUTDIR/step1_basicQC \
  --extract $OUTDIR/step3_prunedSNPs.prune.in \
  --pca 20 \
  --out $OUTDIR/step4_PCA


# ------------------------------------------------------------
# 5. RELATEDNESS (IBD/IBS ESTIMATION)
# ------------------------------------------------------------
echo "[5] Calculating relatedness (IBD/IBS)..."

plink \
  --bfile $OUTDIR/step1_basicQC \
  --genome \
  --out $OUTDIR/step5_relatedness


# ------------------------------------------------------------
# 6. FINAL CLEANED DATASET FOR GWAS
#    (Here no MAF filter yet — apply later depending on your trait)
# ------------------------------------------------------------
echo "[6] Creating final GWAS-ready dataset..."

plink \
  --bfile $OUTDIR/step1_basicQC \
  --make-bed \
  --out $OUTDIR/final_GWAS_ready


echo "GWAS-preparation pipeline complete!"
echo "Final dataset:  $OUTDIR/final_GWAS_ready.{bed,bim,fam}"
echo "PCA results:    $OUTDIR/step4_PCA.eigenvec"
echo "Relatedness:    $OUTDIR/step5_relatedness.genome"
echo "Sex check:      $OUTDIR/step2_sexcheck.sexcheck"
