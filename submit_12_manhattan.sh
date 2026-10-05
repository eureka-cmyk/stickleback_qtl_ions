#!/bin/bash

#SBATCH --job-name=GWAS_Manhattan
#SBATCH --array=1-12%12
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/GWAS_Manhattan_%A_%a.out
#SBATCH --error=/gpfs01/home/shxyt1/slurm_outputs/GWAS_Manhattan_%A_%a.err


set -euo pipefail

module load R-uoneasy/4.4.1-gfbf-2023b

# ============================================================
# 1. R script
#
# Run sbatch from the directory containing:
#
#   GWAS_manhattan_one_trait.R
#   submit_12_manhattan.sh
#
# ============================================================

R_SCRIPT="code/stickleback_qtl_ions/GWAS_manhattan_one_trait.R"


if [ ! -f "${R_SCRIPT}" ]; then

    echo "ERROR: R script not found:"
    echo "${R_SCRIPT}"

    exit 1
fi


# ============================================================
# 2. GWAS input files
#
# Tasks 1-6:
#   5 PCs
#
# Tasks 7-12:
#   no PCs
# ============================================================

MLMA_FILES=(

# ------------------------------------------------------------
# 5 PC
# ------------------------------------------------------------

"/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_nine_main/Ca_no_outlier/Ca_no_outlier_MLMA.mlma"

"/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_nine_main/Co_no_outlier/Co_no_outlier_MLMA.mlma"

"/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_nine_main/Fe_no_outlier/Fe_no_outlier_MLMA.mlma"

"/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_nine_main/Sr_no_outlier/Sr_no_outlier_MLMA.mlma"

"/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_no_outlier/P_log/P_log_MLMA.mlma"

"/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_no_outlier/V_log/V_log_MLMA.mlma"


# ------------------------------------------------------------
# No PC
# ------------------------------------------------------------

"/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_h_w_log_no_outlier/Ca/Ca_MLMA.mlma"

"/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_h_w_log_no_outlier/Co/Co_MLMA.mlma"

"/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_h_w_log_no_outlier/Fe/Fe_MLMA.mlma"

"/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_h_w_log_no_outlier/P/P_MLMA.mlma"

"/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_h_w_log_no_outlier/Sr/Sr_MLMA.mlma"

"/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_h_w_log_no_outlier/V/V_MLMA.mlma"

)


# ============================================================
# 3. Plot titles
# ============================================================

TITLES=(

"Ca_5pc"
"Co_5pc"
"Fe_5pc"
"Sr_5pc"
"P_5pc"
"V_5pc"

"Ca_no_pc"
"Co_no_pc"
"Fe_no_pc"
"P_no_pc"
"Sr_no_pc"
"V_no_pc"

)


# ============================================================
# 4. Convert SLURM task ID to bash array index
#
# SLURM:
#   1-12
#
# Bash:
#   0-11
# ============================================================

IDX=$((SLURM_ARRAY_TASK_ID - 1))


MLMA_FILE="${MLMA_FILES[$IDX]}"

TITLE="${TITLES[$IDX]}"


# ============================================================
# 5. Output directory
# ============================================================

if [ "${SLURM_ARRAY_TASK_ID}" -le 6 ]; then

    OUTPUT_DIR="/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_nine_main/six_main/1.1e-8"

else

    OUTPUT_DIR="/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_h_w_log_no_outlier/1.1e-8"

fi


mkdir -p "${OUTPUT_DIR}"


# ============================================================
# 6. Check input
# ============================================================

if [ ! -f "${MLMA_FILE}" ]; then

    echo "============================================"
    echo "ERROR: MLMA file not found"
    echo "${MLMA_FILE}"
    echo "============================================"

    exit 1

fi


# ============================================================
# 7. Job information
# ============================================================

echo "============================================"
echo "SLURM job ID  : ${SLURM_JOB_ID}"
echo "Array task ID : ${SLURM_ARRAY_TASK_ID}"
echo "Title         : ${TITLE}"
echo "MLMA          : ${MLMA_FILE}"
echo "Output        : ${OUTPUT_DIR}"
echo "Node          : $(hostname)"
echo "Start         : $(date)"
echo "============================================"


# ============================================================
# 8. Run R
# ============================================================

Rscript \
    "${R_SCRIPT}" \
    "${MLMA_FILE}" \
    "${TITLE}" \
    "${OUTPUT_DIR}"


# ============================================================
# 9. Finish
# ============================================================

echo "============================================"
echo "Finished : ${TITLE}"
echo "Time     : $(date)"
echo "============================================"