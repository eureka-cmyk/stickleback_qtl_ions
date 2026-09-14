#!/bin/bash
#SBATCH --job-name=GWAS_manhattan_all
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/GWAS_manhattan_all_%j.log

set -euo pipefail

echo "========================================"
echo "Batch Manhattan and QQ plotting"
echo "Start time: $(date)"
echo "Host: $(hostname)"
echo "SLURM job ID: ${SLURM_JOB_ID}"
echo "========================================"

########################################
# 1. Load R
########################################

module purge
module load R-uoneasy/4.4.1-gfbf-2023b-rstudio

echo "R version:"
R --version | head -1

########################################
# 2. Paths
########################################

PHENO_DIR=/gpfs01/home/shxyt1/data/sticklebacks/pheno/traits/log_no_outlier

GWAS_BASE=/gpfs01/home/shxyt1/data/sticklebacks/gwas/log_no_outlier_22chr_sex_corrected

R_SCRIPT=/gpfs01/home/shxyt1/code/stickleback_qtl_ions/11.3.2V2manhattan.R

WORK_DIR=${GWAS_BASE}/plot_work

COMBINED_DIR=${GWAS_BASE}/combined_plots

TRAIT_LIST=${WORK_DIR}/trait_list.txt

MISSING_LIST=${WORK_DIR}/missing_mlma_files.txt

TRAIT_LOG_DIR=${WORK_DIR}/trait_logs

mkdir -p "${WORK_DIR}"
mkdir -p "${COMBINED_DIR}"
mkdir -p "${TRAIT_LOG_DIR}"

########################################
# 3. Settings
########################################

# Number of R plotting processes running
# simultaneously inside this one SLURM job.
MAX_PARALLEL=4

########################################
# 4. Check R script
########################################

if [[ ! -s "${R_SCRIPT}" ]]; then
    echo "ERROR: R plotting script not found:"
    echo "${R_SCRIPT}"
    exit 1
fi

echo "R script:"
echo "${R_SCRIPT}"

########################################
# 5. Check required R packages
########################################

echo
echo "Checking R packages..."

Rscript -e '
required_packages <- c("data.table", "png")

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  stop(
    "Missing R packages: ",
    paste(missing_packages, collapse = ", ")
  )
}

cat("Required R packages are available.\n")
'

########################################
# 6. Build trait list
########################################

echo
echo "Finding phenotype files in:"
echo "${PHENO_DIR}"

find "${PHENO_DIR}" \
    -maxdepth 1 \
    -type f \
    -name 'pheno_*.txt' \
    -printf '%f\n' \
| sed \
    -e 's/^pheno_//' \
    -e 's/\.txt$//' \
| sort \
> "${TRAIT_LIST}"

N_TRAITS=$(wc -l < "${TRAIT_LIST}")

echo "Phenotypes found: ${N_TRAITS}"
echo

cat -n "${TRAIT_LIST}"

if [[ "${N_TRAITS}" -eq 0 ]]; then
    echo "ERROR: no phenotype files were found."
    exit 1
fi

########################################
# 7. Check all MLMA files before plotting
########################################

echo
echo "Checking MLMA files..."

> "${MISSING_LIST}"

while IFS= read -r TRAIT; do

    MLMA_FILE=${GWAS_BASE}/${TRAIT}/${TRAIT}_MLMA.mlma

    if [[ ! -s "${MLMA_FILE}" ]]; then
        echo "${MLMA_FILE}" >> "${MISSING_LIST}"
    fi

done < "${TRAIT_LIST}"

MISSING_COUNT=$(wc -l < "${MISSING_LIST}")

if [[ "${MISSING_COUNT}" -gt 0 ]]; then

    echo "ERROR: ${MISSING_COUNT} MLMA files are missing or empty:"
    echo

    cat "${MISSING_LIST}"

    echo
    echo "Plotting was not started."
    echo "Missing-file list:"
    echo "${MISSING_LIST}"

    exit 1
fi

echo "All ${N_TRAITS} MLMA files are present."

########################################
# 8. Function for one trait
########################################

run_one_trait() {

    local TRAIT=$1

    local TRAIT_DIR=${GWAS_BASE}/${TRAIT}

    local MLMA_FILE=${TRAIT_DIR}/${TRAIT}_MLMA.mlma

    local OUT_DIR=${TRAIT_DIR}/plots

    local TRAIT_LOG=${TRAIT_LOG_DIR}/${TRAIT}.log

    mkdir -p "${OUT_DIR}"

    {
        echo "========================================"
        echo "Trait: ${TRAIT}"
        echo "Start: $(date)"
        echo "MLMA: ${MLMA_FILE}"
        echo "Output: ${OUT_DIR}"
        echo "========================================"

        Rscript "${R_SCRIPT}" \
            --mode single \
            --base "${GWAS_BASE}" \
            --trait "${TRAIT}" \
            --outdir "${OUT_DIR}"

        MANHATTAN_FILE=${OUT_DIR}/${TRAIT}_Manhattan.png

        QQ_FILE=${OUT_DIR}/${TRAIT}_QQ.png

        if [[ ! -s "${MANHATTAN_FILE}" ]]; then
            echo "ERROR: Manhattan plot missing:"
            echo "${MANHATTAN_FILE}"
            exit 1
        fi

        if [[ ! -s "${QQ_FILE}" ]]; then
            echo "ERROR: QQ plot missing:"
            echo "${QQ_FILE}"
            exit 1
        fi

        echo "Trait completed: ${TRAIT}"
        echo "Finish: $(date)"
        echo "========================================"

    } > "${TRAIT_LOG}" 2>&1
}

########################################
# 9. Run traits in parallel batches
########################################

echo
echo "========================================"
echo "Starting individual plots"
echo "Maximum simultaneous traits: ${MAX_PARALLEL}"
echo "========================================"

declare -a PIDS=()
declare -a PID_TRAITS=()

BATCH_NUMBER=0
COMPLETED_TRAITS=0

wait_for_batch() {

    local FAILED=0

    for INDEX in "${!PIDS[@]}"; do

        PID=${PIDS[$INDEX]}
        TRAIT_NAME=${PID_TRAITS[$INDEX]}

        if wait "${PID}"; then

            COMPLETED_TRAITS=$((COMPLETED_TRAITS + 1))

            echo "Completed ${COMPLETED_TRAITS}/${N_TRAITS}: ${TRAIT_NAME}"

        else

            echo "ERROR: plotting failed for ${TRAIT_NAME}"
            echo "Check log:"
            echo "${TRAIT_LOG_DIR}/${TRAIT_NAME}.log"

            FAILED=1

        fi

    done

    PIDS=()
    PID_TRAITS=()

    if [[ "${FAILED}" -ne 0 ]]; then
        echo "ERROR: at least one plotting task failed."
        exit 1
    fi
}

while IFS= read -r TRAIT; do

    run_one_trait "${TRAIT}" &

    PIDS+=("$!")
    PID_TRAITS+=("${TRAIT}")

    if [[ "${#PIDS[@]}" -ge "${MAX_PARALLEL}" ]]; then

        BATCH_NUMBER=$((BATCH_NUMBER + 1))

        echo
        echo "Waiting for batch ${BATCH_NUMBER}..."

        wait_for_batch

    fi

done < "${TRAIT_LIST}"

# Wait for final incomplete batch
if [[ "${#PIDS[@]}" -gt 0 ]]; then

    BATCH_NUMBER=$((BATCH_NUMBER + 1))

    echo
    echo "Waiting for final batch ${BATCH_NUMBER}..."

    wait_for_batch

fi

echo
echo "All ${COMPLETED_TRAITS} individual plotting tasks completed."

########################################
# 10. Count generated PNG files
########################################

MANHATTAN_COUNT=$(
    find "${GWAS_BASE}" \
        -mindepth 3 \
        -maxdepth 3 \
        -type f \
        -path '*/plots/*_Manhattan.png' \
    | wc -l
)

QQ_COUNT=$(
    find "${GWAS_BASE}" \
        -mindepth 3 \
        -maxdepth 3 \
        -type f \
        -path '*/plots/*_QQ.png' \
    | wc -l
)

echo
echo "Manhattan PNG files found: ${MANHATTAN_COUNT}"
echo "QQ PNG files found: ${QQ_COUNT}"

if [[ "${MANHATTAN_COUNT}" -lt "${N_TRAITS}" ]]; then
    echo "ERROR: fewer Manhattan plots than expected."
    exit 1
fi

if [[ "${QQ_COUNT}" -lt "${N_TRAITS}" ]]; then
    echo "ERROR: fewer QQ plots than expected."
    exit 1
fi

########################################
# 11. Generate combined PDFs
########################################

echo
echo "========================================"
echo "Creating combined PDFs"
echo "========================================"

Rscript "${R_SCRIPT}" \
    --mode combine \
    --base "${GWAS_BASE}" \
    --trait-list "${TRAIT_LIST}" \
    --outdir "${COMBINED_DIR}"

########################################
# 12. Check combined outputs
########################################

ALL_PDF=$(
    find "${COMBINED_DIR}" \
        -maxdepth 1 \
        -type f \
        -name 'All_*_analyses_Manhattan_QQ.pdf' \
        | head -1
)

PREFERRED_PDF=$(
    find "${COMBINED_DIR}" \
        -maxdepth 1 \
        -type f \
        -name 'Preferred_*_elements_Manhattan_QQ.pdf' \
        | head -1
)

SELECTION_FILE=${COMBINED_DIR}/Preferred_analysis_selected_for_each_element.tsv

if [[ -z "${ALL_PDF}" ]] || [[ ! -s "${ALL_PDF}" ]]; then
    echo "ERROR: complete PDF was not generated."
    exit 1
fi

if [[ -z "${PREFERRED_PDF}" ]] || [[ ! -s "${PREFERRED_PDF}" ]]; then
    echo "ERROR: preferred-analysis PDF was not generated."
    exit 1
fi

if [[ ! -s "${SELECTION_FILE}" ]]; then
    echo "ERROR: preferred-analysis selection table was not generated."
    exit 1
fi

########################################
# 13. Final report
########################################

echo
echo "========================================"
echo "All plotting completed successfully"
echo "========================================"
echo "Traits processed: ${N_TRAITS}"
echo "Manhattan plots: ${MANHATTAN_COUNT}"
echo "QQ plots: ${QQ_COUNT}"
echo
echo "All analyses PDF:"
echo "${ALL_PDF}"
echo
echo "Preferred analysis PDF:"
echo "${PREFERRED_PDF}"
echo
echo "Preferred analysis selection table:"
echo "${SELECTION_FILE}"
echo
echo "Individual trait logs:"
echo "${TRAIT_LOG_DIR}"
echo
echo "Finish time: $(date)"
echo "========================================"