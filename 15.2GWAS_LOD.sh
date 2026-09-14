#!/bin/bash
#SBATCH --job-name=GWAS_logLOD
#SBATCH --partition=defq
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --array=1-31
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/gwas_log_LOD_%j_%a.log

echo "Start time: $(date)"
echo "Hostname: $(hostname)"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"

# ============================================================
# 1. Allow overriding via sbatch --export
# ============================================================

GCTA_BIN=${GCTA_BIN:-"/gpfs01/home/shxyt1/app/gcta/gcta-1.95.0-linux-kernel-3-x86_64/squashfs-root/usr/bin/gcta64"}
BFILE=${BFILE:-"/gpfs01/home/shxyt1/data/sticklebacks/PLINK/stickleback_all_clean"}

# NEW: default to log_LOD phenotypes if not provided
PHENO_DIR=${PHENO_DIR:-"/gpfs01/home/shxyt1/data/sticklebacks/pheno/traits/log_LOD"}

# NEW: default to log_LOD PCA covariates (PC1-5) if not provided
QCOVAR=${QCOVAR:-"/gpfs01/home/shxyt1/data/sticklebacks/pheno/traits/log_LOD/log_LOD_PCA_PC1-5.txt"}

# Keep sex as covariate (same as before)
COVAR=${COVAR:-"/gpfs01/home/shxyt1/data/sticklebacks/pheno/traits/sex.txt"}

# NEW: set a new default outroot if not provided
OUTROOT=${OUTROOT:-"/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_log_LOD_5pc_sex"}

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
# 2. Build GRM (only once)
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
# 3. Collect phenotype files (log_LOD)
# IMPORTANT: These files should have NO header line
# ============================================================
PHENO_FILES=($(ls ${PHENO_DIR}/pheno_*_log_LOD.txt))

TOTAL=${#PHENO_FILES[@]}
echo "Total phenotypes: ${TOTAL}"

if [[ ${TOTAL} -lt 1 ]]; then
    echo "Error: no phenotype files found in ${PHENO_DIR}"
    exit 1
fi

if [[ ${SLURM_ARRAY_TASK_ID} -gt ${TOTAL} ]]; then
    echo "Error: array index ${SLURM_ARRAY_TASK_ID} > number of phenotype files (${TOTAL})"
    exit 1
fi

PHENO=${PHENO_FILES[$((SLURM_ARRAY_TASK_ID-1))]}

base=$(basename "${PHENO}")
# Extract element name from: pheno_<EL>_log_LOD.txt
element=$(echo "${base}" | sed 's/^pheno_//; s/_log_LOD\.txt$//')

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
  --reml-no-constrain \
  --reml-maxit 500 \
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
  --reml-no-constrain \
  --reml-maxit 500 \
  --thread-num 8 \
  --out "${OUTDIR}/${element}_MLMA"



echo "Finished ${element} at: $(date)"
