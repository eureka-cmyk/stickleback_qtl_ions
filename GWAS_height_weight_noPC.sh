#!/bin/bash
#SBATCH --job-name=GWAS_HW_noPC
#SBATCH --partition=defq
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --array=1-31
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/GWAS_HW_noPC_%A_%a.log

set -euo pipefail

echo "=================================================="
echo "Start time: $(date)"
echo "Hostname: $(hostname)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "=================================================="

# ============================================================
# 1. Paths
# ============================================================

GCTA_BIN="/gpfs01/home/shxyt1/app/gcta/gcta-1.95.0-linux-kernel-3-x86_64/squashfs-root/usr/bin/gcta64"

BFILE="/gpfs01/home/shxyt1/data/sticklebacks/PLINK/stickleback_all_clean"

PHENO_DIR="/gpfs01/home/shxyt1/data/sticklebacks/pheno/traits/log_no_outlier"

# Quantitative covariates:
# FID IID body_length body_weight
QCOVAR="/gpfs01/home/shxyt1/data/sticklebacks/pheno/traits/height_weight.txt"

# Discrete covariate:
# sex
COVAR="/gpfs01/home/shxyt1/data/sticklebacks/pheno/traits/sex.txt"

# Completely new GWAS output
OUTROOT="/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_h_w_log_no_outlier"

# New GRM for this rerun
GRM_DIR="${OUTROOT}/GRM"
GRM_PREFIX="${GRM_DIR}/stickleback_clean"

mkdir -p "${OUTROOT}"
mkdir -p "${GRM_DIR}"

# ============================================================
# 2. Check input files
# ============================================================

if [[ ! -x "${GCTA_BIN}" ]]; then
    echo "ERROR: GCTA executable not found:"
    echo "${GCTA_BIN}"
    exit 1
fi

if [[ ! -f "${BFILE}.bed" || \
      ! -f "${BFILE}.bim" || \
      ! -f "${BFILE}.fam" ]]; then
    echo "ERROR: PLINK files not found for:"
    echo "${BFILE}"
    exit 1
fi

if [[ ! -f "${QCOVAR}" ]]; then
    echo "ERROR: height/weight covariate file not found:"
    echo "${QCOVAR}"
    exit 1
fi

if [[ ! -f "${COVAR}" ]]; then
    echo "ERROR: sex covariate file not found:"
    echo "${COVAR}"
    exit 1
fi

# ============================================================
# 3. Build a NEW GRM
# ============================================================
#
# Because this is a SLURM array, all 31 jobs start independently.
# We only want ONE task to build the GRM.
#
# mkdir is atomic:
# the first task that successfully creates .grm_build_lock
# becomes the GRM-building task.
# Other tasks wait until the GRM is complete.
# ============================================================

GRM_LOCK="${GRM_DIR}/.grm_build_lock"
GRM_READY="${GRM_DIR}/.grm_ready"
GRM_FAILED="${GRM_DIR}/.grm_failed"

grm_complete() {
    [[ -f "${GRM_PREFIX}.grm.bin" && \
       -f "${GRM_PREFIX}.grm.N.bin" && \
       -f "${GRM_PREFIX}.grm.id" ]]
}

if grm_complete; then

    echo "[GRM] New GRM already exists."

else

    if mkdir "${GRM_LOCK}" 2>/dev/null; then

        echo "[GRM] This array task will build the new GRM."

        rm -f "${GRM_READY}" "${GRM_FAILED}"

        if "${GCTA_BIN}" \
            --bfile "${BFILE}" \
            --make-grm \
            --thread-num 8 \
            --out "${GRM_PREFIX}"
        then

            if grm_complete; then
                touch "${GRM_READY}"
                echo "[GRM] GRM successfully generated."
            else
                echo "[GRM] ERROR: GCTA finished but GRM files are incomplete."
                touch "${GRM_FAILED}"
                rmdir "${GRM_LOCK}" || true
                exit 1
            fi

        else

            echo "[GRM] ERROR: GRM generation failed."
            touch "${GRM_FAILED}"
            rmdir "${GRM_LOCK}" || true
            exit 1

        fi

        rmdir "${GRM_LOCK}" || true

    else

        echo "[GRM] Another array task is building the GRM."
        echo "[GRM] Waiting for completion..."

        # 360 x 30 sec = maximum 3 hours
        for i in $(seq 1 360); do

            if grm_complete; then
                echo "[GRM] GRM is ready."
                break
            fi

            if [[ -f "${GRM_FAILED}" ]]; then
                echo "[GRM] ERROR: GRM generation failed in another task."
                exit 1
            fi

            sleep 30

        done

        if ! grm_complete; then
            echo "[GRM] ERROR: timed out waiting for GRM."
            exit 1
        fi

    fi

fi

# ============================================================
# 4. Determine the 31 elements
# ============================================================
#
# Every element has one pheno_ELEMENT_log.txt.
# *_no_outlier.txt files are replacements for those elements
# where outliers were removed.
#
# Therefore determine the element list ONLY from *_log.txt.
# ============================================================

mapfile -t ELEMENTS < <(
    find "${PHENO_DIR}" \
        -maxdepth 1 \
        -type f \
        -name 'pheno_*_log.txt' \
        -printf '%f\n' \
    | sed 's/^pheno_//; s/_log\.txt$//' \
    | sort
)

TOTAL=${#ELEMENTS[@]}

echo "Total phenotype traits found: ${TOTAL}"

if [[ "${TOTAL}" -ne 31 ]]; then
    echo "ERROR: expected 31 phenotype traits, but found ${TOTAL}."
    echo "Check ${PHENO_DIR}"
    exit 1
fi

if [[ "${SLURM_ARRAY_TASK_ID}" -gt "${TOTAL}" ]]; then
    echo "ERROR: array task ${SLURM_ARRAY_TASK_ID} exceeds ${TOTAL} phenotypes."
    exit 1
fi

element="${ELEMENTS[$((SLURM_ARRAY_TASK_ID - 1))]}"

# ============================================================
# 5. Select phenotype file
# ============================================================
#
# Priority:
#
#   pheno_ELEMENT_no_outlier.txt
#
# if it exists.
#
# Otherwise:
#
#   pheno_ELEMENT_log.txt
#
# ============================================================

NO_OUTLIER_FILE="${PHENO_DIR}/pheno_${element}_no_outlier.txt"
LOG_FILE="${PHENO_DIR}/pheno_${element}_log.txt"

if [[ -f "${NO_OUTLIER_FILE}" ]]; then

    PHENO="${NO_OUTLIER_FILE}"
    PHENO_TYPE="log-transformed + outlier removed"

elif [[ -f "${LOG_FILE}" ]]; then

    PHENO="${LOG_FILE}"
    PHENO_TYPE="log-transformed"

else

    echo "ERROR: phenotype file not found for ${element}"
    exit 1

fi

echo "--------------------------------------------------"
echo "Element        : ${element}"
echo "Phenotype      : ${PHENO}"
echo "Phenotype type : ${PHENO_TYPE}"
echo "QCOVAR         : ${QCOVAR}"
echo "COVAR          : ${COVAR}"
echo "GRM            : ${GRM_PREFIX}"
echo "OUTROOT        : ${OUTROOT}"
echo "--------------------------------------------------"

# ============================================================
# 6. Output directory
# ============================================================

OUTDIR="${OUTROOT}/${element}"
mkdir -p "${OUTDIR}"

# ============================================================
# 7. REML
# ============================================================

echo "[1/2] Running REML for ${element}..."

"${GCTA_BIN}" \
    --grm "${GRM_PREFIX}" \
    --pheno "${PHENO}" \
    --qcovar "${QCOVAR}" \
    --covar "${COVAR}" \
    --reml \
    --thread-num 8 \
    --out "${OUTDIR}/${element}_REML"

echo "[1/2] REML finished for ${element}"

# ============================================================
# 8. MLMA GWAS
# ============================================================

echo "[2/2] Running MLMA GWAS for ${element}..."

"${GCTA_BIN}" \
    --bfile "${BFILE}" \
    --grm "${GRM_PREFIX}" \
    --pheno "${PHENO}" \
    --qcovar "${QCOVAR}" \
    --covar "${COVAR}" \
    --mlma \
    --thread-num 8 \
    --out "${OUTDIR}/${element}_MLMA"

echo "[2/2] MLMA finished for ${element}"

echo "=================================================="
echo "Finished element: ${element}"
echo "End time: $(date)"
echo "=================================================="