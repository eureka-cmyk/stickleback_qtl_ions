#!/bin/bash
#SBATCH --job-name=GWAS_HW_EDA
#SBATCH --partition=defq
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --array=1-31
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/GWAS_noPC_plots/GWAS_HW_EDA_%A_%a.log

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

# Quantitative covariates: height, weight, and EDA class (last column from the EDA file)
HEIGHT_WEIGHT="/gpfs01/home/shxyt1/data/sticklebacks/pheno/traits/height_weight.txt"
EDA_FILE="/gpfs01/home/shxyt1/data/sticklebacks/pheno/traits/pc5_eda_10kb_hetro_class.txt"
QCOVAR="/gpfs01/home/shxyt1/data/sticklebacks/pheno/traits/height_weight_eda.txt"

# Discrete covariate: sex (unchanged)
COVAR="/gpfs01/home/shxyt1/data/sticklebacks/pheno/traits/sex.txt"

# New GWAS output
OUTROOT="/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_h_w_log_eda"

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

if [[ ! -f "${HEIGHT_WEIGHT}" ]]; then
    echo "ERROR: height/weight covariate file not found: ${HEIGHT_WEIGHT}"
    exit 1
fi

if [[ ! -f "${EDA_FILE}" ]]; then
    echo "ERROR: EDA covariate file not found: ${EDA_FILE}"
    exit 1
fi

if [[ ! -f "${COVAR}" ]]; then
    echo "ERROR: sex covariate file not found: ${COVAR}"
    exit 1
fi

# ============================================================
# 3. Create height_weight_eda.txt by matching FID and IID
# ============================================================
# Keep the original height_weight rows and append only the last column
# from pc5_eda_10kb_hetro_class.txt. Matching by IDs avoids relying on
# line order. A lock prevents all array tasks from writing simultaneously.

COVAR_LOCK="${QCOVAR}.${SLURM_ARRAY_JOB_ID}.build_lock"
COVAR_READY="${QCOVAR}.${SLURM_ARRAY_JOB_ID}.ready"
COVAR_FAILED="${QCOVAR}.${SLURM_ARRAY_JOB_ID}.failed"

if [[ ! -f "${COVAR_READY}" ]]; then
    if mkdir "${COVAR_LOCK}" 2>/dev/null; then
        rm -f "${COVAR_READY}" "${COVAR_FAILED}"
        TMP_QCOVAR="${QCOVAR}.tmp.${SLURM_JOB_ID}.${SLURM_ARRAY_TASK_ID}"

        if awk '
            BEGIN { FS="[ \t]+"; OFS="\t"; bad=0 }
            FNR==NR {
                if (NF < 3) { print "ERROR: malformed EDA row at line " FNR > "/dev/stderr"; bad=1; next }
                key=$1 SUBSEP $2
                if (key in eda) { print "ERROR: duplicate FID/IID in EDA file at line " FNR > "/dev/stderr"; bad=1; next }
                eda[key]=$NF
                next
            }
            {
                key=$1 SUBSEP $2
                if (!(key in eda)) { print "ERROR: no EDA value for FID=" $1 ", IID=" $2 > "/dev/stderr"; bad=1; next }
                print $0, eda[key]
                matched++
            }
            END {
                if (matched==0) { print "ERROR: no height/weight rows matched EDA IDs" > "/dev/stderr"; bad=1 }
                if (bad) exit 1
            }
        ' "${EDA_FILE}" "${HEIGHT_WEIGHT}" > "${TMP_QCOVAR}"; then
            mv "${TMP_QCOVAR}" "${QCOVAR}"
            touch "${COVAR_READY}"
            echo "[COVAR] Created ${QCOVAR} by matching FID/IID and appending EDA."
        else
            rm -f "${TMP_QCOVAR}"
            touch "${COVAR_FAILED}"
            rmdir "${COVAR_LOCK}" || true
            exit 1
        fi
        rmdir "${COVAR_LOCK}" || true
    else
        echo "[COVAR] Another array task is creating ${QCOVAR}; waiting..."
        for i in $(seq 1 120); do
            if [[ -f "${COVAR_READY}" && -s "${QCOVAR}" ]]; then
                break
            fi
            if [[ -f "${COVAR_FAILED}" ]]; then
                echo "ERROR: covariate file creation failed in another task."
                exit 1
            fi
            sleep 5
        done
        if [[ ! -f "${COVAR_READY}" || ! -s "${QCOVAR}" ]]; then
            echo "ERROR: timed out waiting for ${QCOVAR}. Check for a stale lock: ${COVAR_LOCK}"
            exit 1
        fi
    fi
else
    echo "[COVAR] Using merged covariate file generated for this array: ${QCOVAR}"
fi

# ============================================================
# 4. Build a NEW GRM
# ============================================================
# Only one array task builds the GRM; the other tasks wait for completion.

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
            --out "${GRM_PREFIX}"; then
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
        echo "[GRM] Another array task is building the GRM; waiting..."
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
# 5. Determine the 31 elements
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
# 6. Select phenotype file
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
# 7. Output directory
# ============================================================

OUTDIR="${OUTROOT}/${element}"
mkdir -p "${OUTDIR}"

# ============================================================
# 8. REML
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
# 9. MLMA GWAS
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
