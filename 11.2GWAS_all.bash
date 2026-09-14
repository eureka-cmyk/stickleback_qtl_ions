#!/bin/bash
#SBATCH --job-name=GWAS
#SBATCH --partition=defq
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --array=1-31
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/gwas_log_no_outlier_%j.log


echo "Start time: $(date)"
echo "Hostname: $(hostname)"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"

# ============================================================
# 1. 接受 sbatch 输入变量（用户可在提交时覆盖）
# ============================================================

GCTA_BIN=${GCTA_BIN:-"/gpfs01/home/shxyt1/app/gcta/gcta-1.95.0-linux-kernel-3-x86_64/squashfs-root/usr/bin/gcta64"}
BFILE=${BFILE:-"/gpfs01/home/shxyt1/data/sticklebacks/PLINK/stickleback_all_clean"}

PHENO_DIR=${PHENO_DIR:? "ERROR: You must provide PHENO_DIR via sbatch --export=PHENO_DIR=..."}
QCOVAR=${QCOVAR:? "ERROR: You must provide QCOVAR via sbatch --export=QCOVAR=..."}
COVAR=${COVAR:-"/gpfs01/home/shxyt1/data/sticklebacks/pheno/traits/sex.txt"}

OUTROOT=${OUTROOT:? "ERROR: You must provide OUTROOT via sbatch --export=OUTROOT=..."}

# GRM will always be in OUTROOT
GRM_DIR="${OUTROOT}/GRM"
mkdir -p "${GRM_DIR}"

echo "--------------------------------------"
echo "PHENO_DIR = ${PHENO_DIR}"
echo "QCOVAR = ${QCOVAR}"
echo "COVAR = ${COVAR}"
echo "OUTROOT = ${OUTROOT}"
echo "GRM_DIR = ${GRM_DIR}"
echo "--------------------------------------"

# ============================================================
# 2. 构建 GRM （只做一次）
# ============================================================
if [[ -f "${GRM_DIR}/stickleback_clean.grm.bin" ]]; then
    echo "[GRM] GRM already exists — skip building."
else
    echo "[GRM] Building GRM..."
    ${GCTA_BIN} \
        --bfile "${BFILE}" \
        --make-grm \
        --thread-num 8 \
        --out "${GRM_DIR}/stickleback_clean"
fi

# ============================================================
# 3. 获取所有 phenotype 文件（自动识别 pheno_*.txt）
# ============================================================
#PHENO_FILES=($(ls ${PHENO_DIR}/pheno_*.txt))
PHENO_FILES=($(ls ${PHENO_DIR}/pheno_*_no_outlier.txt))

TOTAL=${#PHENO_FILES[@]}

echo "Total phenotypes: ${TOTAL}"

if [[ ${SLURM_ARRAY_TASK_ID} -gt ${TOTAL} ]]; then
    echo "Error: array index ${SLURM_ARRAY_TASK_ID} > number of phenotype files (${TOTAL})"
    exit 1
fi

PHENO=${PHENO_FILES[$((SLURM_ARRAY_TASK_ID-1))]}

base=$(basename "${PHENO}")
element=$(echo "${base}" | sed 's/pheno_//; s/.txt//')

echo "Processing phenotype file: ${PHENO}"
echo "Element: ${element}"

OUTDIR="${OUTROOT}/${element}"
mkdir -p "${OUTDIR}"

# ============================================================
# Step 1: REML
# ============================================================
echo "[1] Running REML..."
${GCTA_BIN} \
    --grm "${GRM_DIR}/stickleback_clean" \
    --pheno "${PHENO}" \
    --qcovar "${QCOVAR}" \
    --covar "${COVAR}" \
    --reml \
    --thread-num 8 \
    --out "${OUTDIR}/${element}_REML"

# ============================================================
# Step 2: MLMA GWAS
# ============================================================
echo "[2] Running MLMA GWAS..."
${GCTA_BIN} \
    --bfile "${BFILE}" \
    --grm "${GRM_DIR}/stickleback_clean" \
    --pheno "${PHENO}" \
    --qcovar "${QCOVAR}" \
    --covar "${COVAR}" \
    --mlma \
    --thread-num 8 \
    --out "${OUTDIR}/${element}_MLMA"

echo "Finished ${element} at: $(date)"

#sbatch \--export=PHENO_DIR=/gpfs01/home/shxyt1/data/sticklebacks/pheno/traits/log_no_outlier,\QCOVAR=/gpfs01/home/shxyt1/data/sticklebacks/pheno/traits/height_weight_pc5.txt,\OUTROOT=/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_nine_main \code/stickleback_qtl_ions/11.2GWAS_all.bash 

