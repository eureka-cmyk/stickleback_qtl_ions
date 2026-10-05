#!/bin/bash
#SBATCH --job-name=GWAS_eda_compare
#SBATCH --partition=defq
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/GWAS_eda_compare_%j.log

set -euo pipefail

NO_PC_DIR="/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_h_w_log_eda"
FIVE_PC_DIR="/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_EDA10"
CODE_DIR="/gpfs01/home/shxyt1/code/stickleback_qtl_ions"
R_ONE="${CODE_DIR}/GWAS_compare_one_trait.R"
OUT_DIR="/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_h_w_log_eda/compare"

TRAITS=(V)
mkdir -p "${OUT_DIR}"

if [[ ! -f "${R_ONE}" ]]; then
    echo "ERROR: R script not found: ${R_ONE}"
    exit 1
fi

if [[ ! -d "${NO_PC_DIR}" || ! -d "${FIVE_PC_DIR}" ]]; then
    echo "ERROR: one or both GWAS result directories do not exist."
    echo "No-PC : ${NO_PC_DIR}"
    echo "5-PC  : ${FIVE_PC_DIR}"
    exit 1
fi

# Find one result per element in the 5-PC directory. This accommodates
# phenotype-specific folder names such as Ca_no_outlier and V_log.
find_5pc_mlma() {
    local trait="$1"
    local -a matches=()
    local file dir

    while IFS= read -r file; do
        dir=$(basename "$(dirname "${file}")")
        if [[ "${dir}" =~ ^${trait}(_.*)?$ ]] &&
           [[ "$(basename "${file}")" == "${dir}_MLMA.mlma" ]]; then
            matches+=("${file}")
        fi
    done < <(
        find "${FIVE_PC_DIR}" -mindepth 2 -maxdepth 2 \
            -type f -name '*_MLMA.mlma' | sort
    )

    if [[ "${#matches[@]}" -ne 1 ]]; then
        echo "ERROR: expected exactly one 5-PC MLMA file for ${trait}; found ${#matches[@]}" >&2
        printf '  %s\n' "${matches[@]}" >&2
        return 1
    fi
    printf '%s\n' "${matches[0]}"
}

# Submit the six no-PC and six 5-PC Manhattan plots.
for trait in "${TRAITS[@]}"; do
    no_pc_file="${NO_PC_DIR}/${trait}/${trait}_MLMA.mlma"
    if [[ ! -f "${no_pc_file}" ]]; then
        echo "ERROR: no-PC MLMA file not found: ${no_pc_file}" >&2
        exit 1
    fi

    five_pc_file=$(find_5pc_mlma "${trait}")

    for analysis in 5pc; do
        if [[ "${analysis}" == "no_pc" ]]; then
            mlma_file="${no_pc_file}"
            title="${trait}_eda_no_pc"
        else
            mlma_file="${five_pc_file}"
            title="${trait}_eda_5pc"
        fi

        printf -v RUN_COMMAND \
            'module purge; module load R-uoneasy/4.4.1-gfbf-2023b; Rscript %q %q %q %q' \
            "${R_ONE}" "${mlma_file}" "${title}" "${OUT_DIR}"

        job_id=$(sbatch --parsable \
            --job-name="${title}" \
            --mem=32G \
            --output="/gpfs01/home/shxyt1/slurm_outputs/${title}_%j.log" \
            --wrap="${RUN_COMMAND}")

        echo "Submitted ${title}"
        echo "  MLMA: ${mlma_file}"
        echo "  Job : ${job_id}"
    done
done

echo "All 12 plot jobs submitted. Images will be written to ${OUT_DIR}"
