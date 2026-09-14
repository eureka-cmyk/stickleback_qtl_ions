#!/bin/bash
#SBATCH --job-name=fix_missing_ids_V2
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/fix-missing-ids-V2-%j.log

set -euo pipefail

########################################
# 1. PLINK directory
########################################

PLINK_DIR=/gpfs01/home/shxyt1/data/sticklebacks/PLINK_sex_corrected

cd "${PLINK_DIR}"

shopt -s nullglob

BIM_FILES=(
    stickleback_NC_*_sex_corrected.bim
)

if [[ ${#BIM_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No chromosome BIM files found in:"
    echo "${PLINK_DIR}"
    exit 1
fi

echo "Number of BIM files found: ${#BIM_FILES[@]}"

########################################
# 2. Check and fix missing SNP IDs
########################################

for BIM in "${BIM_FILES[@]}"; do

    echo "========================================"
    echo "Processing: ${BIM}"

    MISSING_BEFORE=$(
        awk '$2=="." || $2=="" {n++} END{print n+0}' "${BIM}"
    )

    echo "Missing SNP IDs before correction: ${MISSING_BEFORE}"

    if [[ "${MISSING_BEFORE}" -gt 0 ]]; then

        TMP_FILE=${BIM}.tmp

        awk '
        BEGIN {
            OFS="\t"
        }

        {
            if ($2=="." || $2=="") {
                $2=$1 ":" $4
            }

            print
        }
        ' "${BIM}" > "${TMP_FILE}"

        if [[ ! -s "${TMP_FILE}" ]]; then
            echo "ERROR: temporary BIM file is empty"
            rm -f "${TMP_FILE}"
            exit 1
        fi

        mv "${TMP_FILE}" "${BIM}"

    else

        echo "No missing SNP IDs. No modification needed."

    fi

    MISSING_AFTER=$(
        awk '$2=="." || $2=="" {n++} END{print n+0}' "${BIM}"
    )

    DUPLICATE_IDS=$(
        awk '{print $2}' "${BIM}" |
        sort |
        uniq -d |
        wc -l
    )

    TOTAL_VARIANTS=$(wc -l < "${BIM}")

    echo "Total variants: ${TOTAL_VARIANTS}"
    echo "Missing SNP IDs after correction: ${MISSING_AFTER}"
    echo "Duplicated SNP IDs: ${DUPLICATE_IDS}"

    if [[ "${MISSING_AFTER}" -ne 0 ]]; then
        echo "ERROR: missing SNP IDs remain in ${BIM}"
        exit 1
    fi

    if [[ "${DUPLICATE_IDS}" -ne 0 ]]; then
        echo "WARNING: duplicated SNP IDs found in ${BIM}"
    fi

done

echo "========================================"
echo "SNP ID checking and correction completed"
echo "PLINK directory: ${PLINK_DIR}"
echo "========================================"