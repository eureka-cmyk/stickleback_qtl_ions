#!/bin/bash
#SBATCH --job-name=MLMA_PCshrink
#SBATCH --partition=defq
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --array=0-5
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/MLMA_PCshrink_%A_%a.log

set -euo pipefail

echo "Start time: $(date)"
echo "Running on: $(hostname)"
echo "Task ID: ${SLURM_ARRAY_TASK_ID}"

# -----------------------------
# Paths
# -----------------------------
GCTA_BIN="/gpfs01/home/shxyt1/app/gcta/gcta-1.95.0-linux-kernel-3-x86_64/squashfs-root/usr/bin/gcta64"
BFILE="/gpfs01/home/shxyt1/data/sticklebacks/PLINK/stickleback_all_clean"

OUTROOT="/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_log_LOD_5pc_sex"
GRM_PREFIX="${OUTROOT}/GRM/stickleback_clean"

PHENO_DIR="/gpfs01/home/shxyt1/data/sticklebacks/pheno/traits/log_LOD"
QCOVAR_5PC="/gpfs01/home/shxyt1/data/sticklebacks/pheno/traits/log_LOD/log_LOD_PCA_PC1-5.txt"
QCOVAR_2PC="/gpfs01/home/shxyt1/data/sticklebacks/pheno/traits/log_LOD/log_LOD_PCA_PC1-2.txt"
COVAR="/gpfs01/home/shxyt1/data/sticklebacks/pheno/traits/sex.txt"

# Elements and models
elements=(Ca Sr Ti)
models=("sexOnly" "PC12")   # 0=sexOnly, 1=PC12

# -----------------------------
# Check GRM exists (do not rebuild)
# -----------------------------
for f in "${GRM_PREFIX}.grm.bin" "${GRM_PREFIX}.grm.id" "${GRM_PREFIX}.grm.N.bin"; do
  if [[ ! -f "$f" ]]; then
    echo "[ABORT] Missing GRM file: $f"
    exit 1
  fi
done

# -----------------------------
# Build PC1-2 file once (race-safe)
# Only task 0 creates it; others wait a bit if needed
# -----------------------------
if [[ ${SLURM_ARRAY_TASK_ID} -eq 0 ]]; then
  if [[ ! -f "$QCOVAR_2PC" ]]; then
    echo "[INFO] Creating PC1-2 covariate file: $QCOVAR_2PC"
    awk 'BEGIN{OFS="\t"} NR==1{print $1,$2,$3,$4; next} {print $1,$2,$3,$4}' \
      "$QCOVAR_5PC" > "$QCOVAR_2PC"
  fi
else
  # wait up to ~60s for task0 to create it
  for i in $(seq 1 30); do
    [[ -f "$QCOVAR_2PC" ]] && break
    sleep 2
  done
fi

# -----------------------------
# Map task ID -> (element, model)
# task 0-5:
#   0: Ca sexOnly
#   1: Ca PC12
#   2: Sr sexOnly
#   3: Sr PC12
#   4: Ti sexOnly
#   5: Ti PC12
# -----------------------------
tid=${SLURM_ARRAY_TASK_ID}
el_idx=$(( tid / 2 ))
m_idx=$(( tid % 2 ))

el="${elements[$el_idx]}"
model="${models[$m_idx]}"

PHENO="${PHENO_DIR}/pheno_${el}_log_LOD.txt"
OUTDIR="${OUTROOT}/${el}"
mkdir -p "$OUTDIR"

if [[ ! -f "$PHENO" ]]; then
  echo "[ABORT] Missing phenotype file: $PHENO"
  exit 1
fi

echo "Element: $el"
echo "Model:   $model"
echo "PHENO:   $PHENO"

# -----------------------------
# Run GCTA MLMA
# -----------------------------
if [[ "$model" == "sexOnly" ]]; then
  OUT="${OUTDIR}/${el}_MLMA_sexOnly"
  ${GCTA_BIN} \
    --bfile "${BFILE}" \
    --grm "${GRM_PREFIX}" \
    --pheno "${PHENO}" \
    --covar "${COVAR}" \
    --mlma \
    --reml-no-constrain \
    --reml-maxit 2000 \
    --thread-num 8 \
    --out "${OUT}"

else
  OUT="${OUTDIR}/${el}_MLMA_PC12"
  ${GCTA_BIN} \
    --bfile "${BFILE}" \
    --grm "${GRM_PREFIX}" \
    --pheno "${PHENO}" \
    --qcovar "${QCOVAR_2PC}" \
    --covar "${COVAR}" \
    --mlma \
    --reml-no-constrain \
    --reml-maxit 2000 \
    --thread-num 8 \
    --out "${OUT}"
fi

echo "Done: $(date)"
echo "Output prefix: ${OUT}"
