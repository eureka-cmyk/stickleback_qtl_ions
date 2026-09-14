#!/bin/bash
#SBATCH --job-name=plink_merge_V2
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/plink-merge-V2-%j.out

set -euo pipefail

module load plink-uoneasy/1.9b_6.21-x86_64

########################################
# 1. Paths
########################################

GENOME_NAME=GCF_016920845.1_GAculeatus_UGA_version5_genomic

CHROMOSOME_BED=/share/MacColl_shared/QTL/genomes/${GENOME_NAME}.bed

PLINK_DIR=/gpfs01/home/shxyt1/data/sticklebacks/PLINK_sex_corrected

OUTPUT_PREFIX=${PLINK_DIR}/stickleback_allchr_sex_corrected

PREFIX_LIST=${PLINK_DIR}/all_plink_prefixes_sex_corrected.txt

MERGE_LIST=${PLINK_DIR}/merge_list_sex_corrected.txt

cd "${PLINK_DIR}"

echo "Working directory:"
echo "$(pwd)"

########################################
# 2. Build chromosome list in genome order
########################################

> "${PREFIX_LIST}"

while read -r chr rest; do

    prefix=${PLINK_DIR}/stickleback_${chr}_sex_corrected

    for ext in bed bim fam; do
        if [[ ! -s "${prefix}.${ext}" ]]; then
            echo "ERROR: missing or empty PLINK file:"
            echo "${prefix}.${ext}"
            exit 1
        fi
    done

    echo "${prefix}" >> "${PREFIX_LIST}"

done < <(head -22 "${CHROMOSOME_BED}")

NUMBER_PREFIXES=$(wc -l < "${PREFIX_LIST}")

echo "Number of chromosome datasets found: ${NUMBER_PREFIXES}"

if [[ "${NUMBER_PREFIXES}" -ne 22 ]]; then
    echo "ERROR: expected 22 chromosome datasets"
    exit 1
fi

echo "Chromosome datasets will be merged in this order:"
cat "${PREFIX_LIST}"

########################################
# 3. Check FAM files
########################################

BASE_PREFIX=$(head -1 "${PREFIX_LIST}")

while read -r prefix; do

    if ! cmp -s "${BASE_PREFIX}.fam" "${prefix}.fam"; then
        echo "ERROR: FAM file differs from the base dataset:"
        echo "${prefix}.fam"
        exit 1
    fi

done < "${PREFIX_LIST}"

echo "FAM file check: PASS"
echo "All chromosome files contain identical samples in identical order."

########################################
# 4. Prepare merge list
########################################

tail -n +2 "${PREFIX_LIST}" > "${MERGE_LIST}"

echo "Base dataset:"
echo "${BASE_PREFIX}"

echo "Merge list:"
cat "${MERGE_LIST}"

########################################
# 5. Prevent accidental overwriting
########################################

for ext in bed bim fam log nosex; do
    if [[ -e "${OUTPUT_PREFIX}.${ext}" ]]; then
        echo "ERROR: output already exists:"
        echo "${OUTPUT_PREFIX}.${ext}"
        echo "Move or delete the existing output before rerunning."
        exit 1
    fi
done

########################################
# 6. Merge 22 chromosome datasets
########################################

echo "Starting PLINK merge..."

plink \
    --bfile "${BASE_PREFIX}" \
    --merge-list "${MERGE_LIST}" \
    --allow-extra-chr \
    --keep-allele-order \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --memory 30000 \
    --make-bed \
    --out "${OUTPUT_PREFIX}"

########################################
# 7. Check merged output
########################################

for ext in bed bim fam; do

    if [[ ! -s "${OUTPUT_PREFIX}.${ext}" ]]; then
        echo "ERROR: merged output is missing or empty:"
        echo "${OUTPUT_PREFIX}.${ext}"
        exit 1
    fi

done

SAMPLE_COUNT=$(wc -l < "${OUTPUT_PREFIX}.fam")
VARIANT_COUNT=$(wc -l < "${OUTPUT_PREFIX}.bim")

MISSING_IDS=$(
    awk '$2=="." || $2=="" {n++} END{print n+0}' \
    "${OUTPUT_PREFIX}.bim"
)

DUPLICATE_IDS=$(
    awk '{print $2}' "${OUTPUT_PREFIX}.bim" |
    sort |
    uniq -d |
    wc -l
)

########################################
# 8. Final report
########################################

echo "========================================"
echo "PLINK merge completed successfully"
echo "Samples: ${SAMPLE_COUNT}"
echo "Variants: ${VARIANT_COUNT}"
echo "Missing variant IDs: ${MISSING_IDS}"
echo "Duplicated variant IDs: ${DUPLICATE_IDS}"
echo
echo "Merged files:"
echo "${OUTPUT_PREFIX}.bed"
echo "${OUTPUT_PREFIX}.bim"
echo "${OUTPUT_PREFIX}.fam"
echo "========================================"