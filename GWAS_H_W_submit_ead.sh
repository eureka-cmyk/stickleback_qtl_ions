#!/bin/bash
set -euo pipefail

# EDA-adjusted GWAS result directory
BASE_DIR="/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_h_w_log_eda"
CODE_DIR="/gpfs01/home/shxyt1/code/stickleback_qtl_ions"
R_ONE="${CODE_DIR}/GWAS_H_W_MAH_one_trait_eda.R"
LOG_DIR="/gpfs01/home/shxyt1/slurm_outputs/GWAS_noPC_plots"
PLOT_DIR="${BASE_DIR}/plots_two_thresholds"

mkdir -p "${LOG_DIR}" "${PLOT_DIR}"

if [[ ! -f "${R_ONE}" ]]; then
    echo "ERROR: R script not found: ${R_ONE}"
    exit 1
fi

# Submit one plotting job for each existing MLMA result.
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
echo "EDA-adjusted MLMA files found: ${N}"
echo "============================================"

if [[ "${N}" -eq 0 ]]; then
    echo "ERROR: no MLMA files found under ${BASE_DIR}"
    exit 1
fi

JOB_IDS=()
for MLMA in "${MLMA_FILES[@]}"; do
    FILE_NAME=$(basename "${MLMA}")
    TRAIT=${FILE_NAME%_MLMA.mlma}

    printf -v RUN_COMMAND \
        'module purge; module load R-uoneasy/4.4.1-gfbf-2023b; Rscript %q %q %q %q' \
        "${R_ONE}" \
        "${MLMA}" \
        "${TRAIT}" \
        "${PLOT_DIR}"

    JOB_ID=$(
        sbatch \
            --parsable \
            --job-name="EDAplot_${TRAIT}" \
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

echo ""
echo "Submitted ${#JOB_IDS[@]} plotting jobs."
echo "Job IDs:"
printf '%s\n' "${JOB_IDS[@]}"
echo "Plot output directory: ${PLOT_DIR}"
echo "Each successful trait produces a Manhattan PNG and a QQ PNG."
