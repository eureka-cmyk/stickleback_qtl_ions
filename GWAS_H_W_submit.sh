#!/bin/bash

set -euo pipefail


# ============================================================
# Paths
# ============================================================

BASE_DIR="/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_h_w_log_no_outlier"

CODE_DIR="/gpfs01/home/shxyt1/code/stickleback_qtl_ions"

R_ONE="${CODE_DIR}/GWAS_H_W_MAH_one_trait.R"

LOG_DIR="/gpfs01/home/shxyt1/slurm_outputs/GWAS_noPC_plots"


mkdir -p "${LOG_DIR}"

mkdir -p "${BASE_DIR}/plots_two_thresholds"

mkdir -p "${BASE_DIR}/per_trait_tables"


# ============================================================
# Check R script
# ============================================================

if [[ ! -f "${R_ONE}" ]]; then

    echo "ERROR: R script not found:"
    echo "${R_ONE}"

    exit 1

fi


# ============================================================
# Find all successful MLMA files
#
# Only existing GWAS results are submitted.
# Traits whose GWAS failed are automatically ignored.
# ============================================================

mapfile -t MLMA_FILES < <(

    find "${BASE_DIR}" \
        -mindepth 2 \
        -maxdepth 2 \
        -type f \
        -name '*_MLMA.mlma' \
        | sort

)


N=${#MLMA_FILES[@]}


echo "============================================"
echo "Successful MLMA files found: ${N}"
echo "============================================"


if [[ "${N}" -eq 0 ]]; then

    echo "ERROR: no MLMA files found."

    exit 1

fi


# ============================================================
# Submit one independent plotting job per trait
# ============================================================

JOB_IDS=()


for MLMA in "${MLMA_FILES[@]}"; do

    FILE_NAME=$(basename "${MLMA}")

    TRAIT=${FILE_NAME%_MLMA.mlma}


    echo ""
    echo "Submitting trait: ${TRAIT}"
    echo "MLMA: ${MLMA}"


    # --------------------------------------------------------
    # Command executed inside each SLURM job
    # --------------------------------------------------------

    printf -v RUN_COMMAND \
        'module purge; module load R-uoneasy/4.4.1-gfbf-2023b; Rscript %q %q %q' \
        "${R_ONE}" \
        "${MLMA}" \
        "${BASE_DIR}"


    # --------------------------------------------------------
    # Submit
    # --------------------------------------------------------

    JOB_ID=$(

        sbatch \
            --parsable \
            --job-name="MAH_${TRAIT}" \
            --partition=defq \
            --cpus-per-task=1 \
            --mem=12G \
            --time=02:00:00 \
            --output="${LOG_DIR}/${TRAIT}_plot_%j.out" \
            --wrap="${RUN_COMMAND}"

    )


    JOB_IDS+=("${JOB_ID}")


    echo "Submitted ${TRAIT}: job ${JOB_ID}"

done


# ============================================================
# Finished submission
# ============================================================

echo ""
echo "============================================"
echo "Submitted ${#JOB_IDS[@]} plotting jobs."
echo "============================================"

echo ""
echo "Job IDs:"

printf '%s\n' "${JOB_IDS[@]}"

echo ""
echo "Plot output directory:"
echo "${BASE_DIR}/plots_two_thresholds"

echo ""
echo "Per-trait table directory:"
echo "${BASE_DIR}/per_trait_tables"

echo ""
echo "After these jobs finish, run:"
echo ""
echo "sbatch ${CODE_DIR}/GWAS_noPC_merge_manhattan.sh"

echo ""
echo "============================================"
echo "Plot submission complete."
echo "============================================"