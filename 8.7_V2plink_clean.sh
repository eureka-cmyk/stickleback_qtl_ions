#!/bin/bash
#SBATCH --job-name=plink_clean_V2
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/plink-clean-V2-%j.log

set -euo pipefail

module load plink-uoneasy/1.9b_6.21-x86_64


########################################
# 1. Input and output paths
########################################

PLINK_DIR=/gpfs01/home/shxyt1/data/sticklebacks/PLINK_sex_corrected

PLINK_MASTER=${PLINK_DIR}/stickleback_allchr_sex_corrected

# Two-column file:
# QTL22002    F
# QTL22005    M
GSEX_FILE=/gpfs01/home/shxyt1/data/sticklebacks/VCF/ploidy/Gsex_bcftools.tsv

# Per-chromosome Y dataset; used to obtain Y variant IDs
Y_PREFIX=${PLINK_DIR}/stickleback_NC_053233.1_sex_corrected

# Separate output directory
OUT_DIR=/gpfs01/home/shxyt1/data/sticklebacks/PLINK_clean_sex_corrected

mkdir -p "${OUT_DIR}"


########################################
# 2. Output prefixes
########################################

SEX_UPDATED=${OUT_DIR}/stickleback_allchr_sex_updated

MAIN_CLEAN=${OUT_DIR}/stickleback_all_clean_sex_corrected

Y_CLEAN=${OUT_DIR}/stickleback_Y_clean_sex_corrected

PLINK_SEX_FILE=${OUT_DIR}/sample_sex_for_plink.txt

MALE_KEEP_FILE=${OUT_DIR}/male_samples_for_plink.txt

Y_VARIANT_LIST=${OUT_DIR}/Y_variant_ids.txt


########################################
# 3. Check required input files
########################################

for EXT in bed bim fam; do

    if [[ ! -s "${PLINK_MASTER}.${EXT}" ]]; then
        echo "ERROR: missing input file:"
        echo "${PLINK_MASTER}.${EXT}"
        exit 1
    fi

    if [[ ! -s "${Y_PREFIX}.${EXT}" ]]; then
        echo "ERROR: missing Y chromosome file:"
        echo "${Y_PREFIX}.${EXT}"
        exit 1
    fi

done

if [[ ! -s "${GSEX_FILE}" ]]; then
    echo "ERROR: missing sex file:"
    echo "${GSEX_FILE}"
    exit 1
fi


########################################
# 4. Create PLINK sex file
#
# PLINK format:
# FID IID SEX
#
# SEX:
# 1 = male
# 2 = female
########################################

awk '
BEGIN {
    OFS="\t"
}

$2=="M" {
    print $1, $1, 1
}

$2=="F" {
    print $1, $1, 2
}
' "${GSEX_FILE}" > "${PLINK_SEX_FILE}"

SEX_TOTAL=$(wc -l < "${PLINK_SEX_FILE}")

SEX_MALES=$(
    awk '$3==1 {n++} END {print n+0}' \
    "${PLINK_SEX_FILE}"
)

SEX_FEMALES=$(
    awk '$3==2 {n++} END {print n+0}' \
    "${PLINK_SEX_FILE}"
)

echo "Samples in PLINK sex file: ${SEX_TOTAL}"
echo "Males: ${SEX_MALES}"
echo "Females: ${SEX_FEMALES}"

if [[ "${SEX_TOTAL}" -ne 352 ]]; then
    echo "ERROR: expected 352 samples"
    exit 1
fi

if [[ "${SEX_MALES}" -ne 103 ]]; then
    echo "ERROR: expected 103 males"
    exit 1
fi

if [[ "${SEX_FEMALES}" -ne 249 ]]; then
    echo "ERROR: expected 249 females"
    exit 1
fi


########################################
# 5. Create male-only keep file
########################################

awk '
BEGIN {
    OFS="\t"
}

$3==1 {
    print $1, $2
}
' "${PLINK_SEX_FILE}" > "${MALE_KEEP_FILE}"

if [[ $(wc -l < "${MALE_KEEP_FILE}") -ne 103 ]]; then
    echo "ERROR: male keep file should contain 103 samples"
    exit 1
fi


########################################
# 6. Add sex information to merged PLINK
########################################

echo "========================================"
echo "[1] Updating sex information"
echo "========================================"

plink \
    --bfile "${PLINK_MASTER}" \
    --update-sex "${PLINK_SEX_FILE}" \
    --allow-extra-chr \
    --keep-allele-order \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --memory 30000 \
    --make-bed \
    --out "${SEX_UPDATED}"


########################################
# 7. Check updated FAM
########################################

UPDATED_TOTAL=$(wc -l < "${SEX_UPDATED}.fam")

UPDATED_MALES=$(
    awk '$5==1 {n++} END {print n+0}' \
    "${SEX_UPDATED}.fam"
)

UPDATED_FEMALES=$(
    awk '$5==2 {n++} END {print n+0}' \
    "${SEX_UPDATED}.fam"
)

UPDATED_UNKNOWN=$(
    awk '$5==0 {n++} END {print n+0}' \
    "${SEX_UPDATED}.fam"
)

echo "Updated FAM samples: ${UPDATED_TOTAL}"
echo "Updated FAM males: ${UPDATED_MALES}"
echo "Updated FAM females: ${UPDATED_FEMALES}"
echo "Updated FAM unknown sex: ${UPDATED_UNKNOWN}"

if [[ "${UPDATED_TOTAL}" -ne 352 ||
      "${UPDATED_MALES}" -ne 103 ||
      "${UPDATED_FEMALES}" -ne 249 ||
      "${UPDATED_UNKNOWN}" -ne 0 ]]; then

    echo "ERROR: sex information was not updated correctly"
    exit 1
fi


########################################
# 8. Create Y variant list
########################################

awk '{print $2}' \
    "${Y_PREFIX}.bim" \
| sort -u \
> "${Y_VARIANT_LIST}"

Y_VARIANT_COUNT=$(wc -l < "${Y_VARIANT_LIST}")

echo "Y chromosome variants identified: ${Y_VARIANT_COUNT}"

if [[ "${Y_VARIANT_COUNT}" -eq 0 ]]; then
    echo "ERROR: Y variant list is empty"
    exit 1
fi


########################################
# 9. Main cleaning
#
# Exclude Y chromosome first.
# Retain autosomes and X/chr19.
#
# --geno 0.05:
# remove variants with >5% missingness
#
# --maf 0.01:
# remove variants with MAF <1%
########################################

echo "========================================"
echo "[2] Cleaning autosomes and X chromosome"
echo "========================================"

plink \
    --bfile "${SEX_UPDATED}" \
    --exclude "${Y_VARIANT_LIST}" \
    --geno 0.05 \
    --maf 0.01 \
    --allow-extra-chr \
    --keep-allele-order \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --memory 30000 \
    --make-bed \
    --out "${MAIN_CLEAN}"


########################################
# 10. Clean Y chromosome separately
#
# Only male samples are included.
########################################

echo "========================================"
echo "[3] Cleaning Y chromosome in males only"
echo "========================================"

plink \
    --bfile "${SEX_UPDATED}" \
    --extract "${Y_VARIANT_LIST}" \
    --keep "${MALE_KEEP_FILE}" \
    --geno 0.05 \
    --maf 0.01 \
    --allow-extra-chr \
    --keep-allele-order \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --memory 30000 \
    --make-bed \
    --out "${Y_CLEAN}"


########################################
# 11. Check outputs
########################################

for PREFIX in "${MAIN_CLEAN}" "${Y_CLEAN}"; do

    for EXT in bed bim fam; do

        if [[ ! -s "${PREFIX}.${EXT}" ]]; then
            echo "ERROR: missing or empty output:"
            echo "${PREFIX}.${EXT}"
            exit 1
        fi

    done

done


########################################
# 12. Summaries
########################################

MAIN_SAMPLES=$(wc -l < "${MAIN_CLEAN}.fam")
MAIN_VARIANTS=$(wc -l < "${MAIN_CLEAN}.bim")

Y_SAMPLES=$(wc -l < "${Y_CLEAN}.fam")
Y_VARIANTS=$(wc -l < "${Y_CLEAN}.bim")

MAIN_MISSING_IDS=$(
    awk '$2=="." || $2=="" {n++} END {print n+0}' \
    "${MAIN_CLEAN}.bim"
)

MAIN_DUPLICATE_IDS=$(
    awk '{print $2}' "${MAIN_CLEAN}.bim" |
    sort |
    uniq -d |
    wc -l
)


########################################
# 13. Final report
########################################

echo "========================================"
echo "PLINK cleaning completed successfully"
echo "========================================"

echo "Main autosome + X dataset:"
echo "Samples: ${MAIN_SAMPLES}"
echo "Variants: ${MAIN_VARIANTS}"
echo "Missing variant IDs: ${MAIN_MISSING_IDS}"
echo "Duplicated variant IDs: ${MAIN_DUPLICATE_IDS}"
echo "Files:"
echo "${MAIN_CLEAN}.{bed,bim,fam}"

echo

echo "Male-only Y dataset:"
echo "Samples: ${Y_SAMPLES}"
echo "Variants: ${Y_VARIANTS}"
echo "Files:"
echo "${Y_CLEAN}.{bed,bim,fam}"

echo "========================================"