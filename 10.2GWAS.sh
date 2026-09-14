#!/bin/bash
#SBATCH --job-name=Na_RINT_GWAS
#SBATCH --partition=defq
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/Na_RINT_GWAS_%j.log

echo "Start time: $(date)"
echo "Running on: $(hostname)"

# === 1. Paths and basic settings ============================================
GCTA_BIN="/gpfs01/home/shxyt1/app/gcta/gcta-1.95.0-linux-kernel-3-x86_64/squashfs-root/usr/bin/gcta64"

BFILE="/gpfs01/home/shxyt1/data/sticklebacks/PLINK/stickleback_all_clean"

PHENO="/gpfs01/home/shxyt1/data/sticklebacks/pheno/traits/pheno_Na_RINT.txt"
QCOVAR="/gpfs01/home/shxyt1/data/sticklebacks/pheno/traits/pcs5.txt"
COVAR="/gpfs01/home/shxyt1/data/sticklebacks/pheno/traits/sex.txt"

OUTDIR="/gpfs01/home/shxyt1/data/sticklebacks/gwas/Na_RINT"
mkdir -p "${OUTDIR}"

echo "Genotype  : ${BFILE}"
echo "Phenotype : ${PHENO}"
echo "Qcovar    : ${QCOVAR}"
echo "Covar     : ${COVAR}"
echo "Outdir    : ${OUTDIR}"

# === 2. Step 1: Make GRM ====================================================
# GRM 用 stickleback_all_clean 全部 SNP
echo "[1] Making GRM from cleaned PLINK data..."
${GCTA_BIN} \
  --bfile "${BFILE}" \
  --make-grm \
  --thread-num 8 \
  --out "${OUTDIR}/stickleback_all_clean"

echo "[1] GRM files:"
echo "  ${OUTDIR}/stickleback_all_clean.grm.bin"
echo "  ${OUTDIR}/stickleback_all_clean.grm.N.bin"
echo "  ${OUTDIR}/stickleback_all_clean.grm.id"

# === 3. Step 2 (optional but nice): REML to estimate h2 =====================
# 这一步是看一下 Na_RINT 的遗传率，大概有多少 signal
# 不想跑可以直接注释掉整个 block
echo "[2] Running REML for Na_RINT (heritability estimation)..."
${GCTA_BIN} \
  --grm "${OUTDIR}/stickleback_all_clean" \
  --pheno "${PHENO}" \
  --qcovar "${QCOVAR}" \
  --covar "${COVAR}" \
  --reml \
  --thread-num 8 \
  --out "${OUTDIR}/Na_RINT_REML"

echo "[2] REML done. Check:"
echo "  ${OUTDIR}/Na_RINT_REML.log"
echo "  ${OUTDIR}/Na_RINT_REML.hsq"

# === 4. Step 3: MLMA GWAS ===================================================
# 线性混合模型关联分析 (MLMA)，使用 GRM 建模个体间 relatedness
echo "[3] Running MLMA GWAS for Na_RINT..."
${GCTA_BIN} \
  --bfile "${BFILE}" \
  --grm "${OUTDIR}/stickleback_all_clean" \
  --pheno "${PHENO}" \
  --qcovar "${QCOVAR}" \
  --covar "${COVAR}" \
  --mlma \
  --thread-num 8 \
  --out "${OUTDIR}/Na_RINT_MLMA"

echo "[3] MLMA done. Main results:"
echo "  ${OUTDIR}/Na_RINT_MLMA.mlma"

echo "Finished at: $(date)"
