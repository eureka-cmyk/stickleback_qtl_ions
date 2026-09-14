#!/bin/bash
#SBATCH --job-name=HeightWeight_GWAS
#SBATCH --partition=defq
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/HeightWeight_GWAS_%j.log

echo "Start time: $(date)"
echo "Running on: $(hostname)"

# === 1. Paths and settings ==================================================
GCTA_BIN="/gpfs01/home/shxyt1/app/gcta/gcta-1.95.0-linux-kernel-3-x86_64/squashfs-root/usr/bin/gcta64"

BFILE="/gpfs01/home/shxyt1/data/sticklebacks/PLINK/stickleback_all_clean"

PHENO_HEIGHT="/gpfs01/home/shxyt1/data/sticklebacks/pheno/traits/pheno_height.txt"
PHENO_WEIGHT="/gpfs01/home/shxyt1/data/sticklebacks/pheno/traits/pheno_weight.txt"
QCOVAR="/gpfs01/home/shxyt1/data/sticklebacks/pheno/traits/pcs5.txt"
COVAR="/gpfs01/home/shxyt1/data/sticklebacks/pheno/traits/sex.txt"

OUTDIR_BASE="/gpfs01/home/shxyt1/data/sticklebacks/gwas/height_weight"
mkdir -p "${OUTDIR_BASE}"

echo "Genotype file : ${BFILE}"
echo "Height pheno  : ${PHENO_HEIGHT}"
echo "Weight pheno  : ${PHENO_WEIGHT}"
echo "Qcovar file   : ${QCOVAR}"
echo "Covar file    : ${COVAR}"
echo "Output base   : ${OUTDIR_BASE}"

# === 2. Make GRM once =======================================================
echo "[1] Making GRM from cleaned PLINK data..."
${GCTA_BIN} \
  --bfile "${BFILE}" \
  --make-grm \
  --thread-num 8 \
  --out "${OUTDIR_BASE}/stickleback_all_clean"

echo "[1] GRM finished."

# === 3. Height GWAS =========================================================
HEIGHT_OUTDIR="${OUTDIR_BASE}/height"
mkdir -p "${HEIGHT_OUTDIR}"

echo "[2] Running REML for height..."
${GCTA_BIN} \
  --grm "${OUTDIR_BASE}/stickleback_all_clean" \
  --pheno "${PHENO_HEIGHT}" \
  --qcovar "${QCOVAR}" \
  --covar "${COVAR}" \
  --reml \
  --thread-num 8 \
  --out "${HEIGHT_OUTDIR}/height_REML"

echo "[3] Running MLMA GWAS for height..."
${GCTA_BIN} \
  --bfile "${BFILE}" \
  --grm "${OUTDIR_BASE}/stickleback_all_clean" \
  --pheno "${PHENO_HEIGHT}" \
  --qcovar "${QCOVAR}" \
  --covar "${COVAR}" \
  --mlma \
  --thread-num 8 \
  --out "${HEIGHT_OUTDIR}/height_MLMA"

echo "[3] Height GWAS finished."

# === 4. Weight GWAS =========================================================
WEIGHT_OUTDIR="${OUTDIR_BASE}/weight"
mkdir -p "${WEIGHT_OUTDIR}"

echo "[4] Running REML for weight..."
${GCTA_BIN} \
  --grm "${OUTDIR_BASE}/stickleback_all_clean" \
  --pheno "${PHENO_WEIGHT}" \
  --qcovar "${QCOVAR}" \
  --covar "${COVAR}" \
  --reml \
  --thread-num 8 \
  --out "${WEIGHT_OUTDIR}/weight_REML"

echo "[5] Running MLMA GWAS for weight..."
${GCTA_BIN} \
  --bfile "${BFILE}" \
  --grm "${OUTDIR_BASE}/stickleback_all_clean" \
  --pheno "${PHENO_WEIGHT}" \
  --qcovar "${QCOVAR}" \
  --covar "${COVAR}" \
  --mlma \
  --thread-num 8 \
  --out "${WEIGHT_OUTDIR}/weight_MLMA"

echo "[5] Weight GWAS finished."

echo "Finished at: $(date)"