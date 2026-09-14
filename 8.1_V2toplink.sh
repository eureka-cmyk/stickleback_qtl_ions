#!/bin/bash
#SBATCH --job-name=plink_bcf2bed_V2
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --array=1-22%4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/plink-bcf2bed-V2-%A_%a.out

set -euo pipefail

module load bcftools-uoneasy/1.19-GCC-13.2.0
module load plink-uoneasy/1.9b_6.21-x86_64

########################################
# Paths
########################################

GENOME_NAME=GCF_016920845.1_GAculeatus_UGA_version5_genomic

CHROMOSOME_BED=/share/MacColl_shared/QTL/genomes/${GENOME_NAME}.bed

BCF_DIR=/gpfs01/home/shxyt1/data/sticklebacks/VCF/ploidy_aware_HWEPops_MQ10_BQ20_sex_corrected

CHR_MAP=/gpfs01/home/shxyt1/data/sticklebacks/chr_map.txt

OUT_DIR=/gpfs01/home/shxyt1/data/sticklebacks/PLINK_sex_corrected

TEMP_BASE=/gpfs01/home/shxyt1/temp/plink_bcf2bed_V2

mkdir -p "${OUT_DIR}"
mkdir -p "${TEMP_BASE}"

########################################
# Select chromosome
########################################

CHR=$(awk -v task="${SLURM_ARRAY_TASK_ID}" '
NR == task {
    print $1
    exit
}
' "${CHROMOSOME_BED}")

if [[ -z "${CHR}" ]]; then
    echo "ERROR: chromosome not found for task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

INPUT_BCF=${BCF_DIR}/stickleback_${CHR}_sorted.bcf

OUTPUT_PREFIX=${OUT_DIR}/stickleback_${CHR}_sex_corrected

TEMP_DIR=${TEMP_BASE}/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

TEMP_BCF=${TEMP_DIR}/stickleback_${CHR}_biallelic_snps_renamed.bcf

mkdir -p "${TEMP_DIR}"

trap 'rm -rf "${TEMP_DIR}"' EXIT

########################################
# Input checks
########################################

if [[ ! -s "${INPUT_BCF}" ]]; then
    echo "ERROR: input BCF missing or empty:"
    echo "${INPUT_BCF}"
    exit 1
fi

if [[ ! -s "${INPUT_BCF}.csi" ]]; then
    echo "ERROR: input BCF index missing:"
    echo "${INPUT_BCF}.csi"
    exit 1
fi

if [[ ! -s "${CHR_MAP}" ]]; then
    echo "ERROR: chromosome map missing or empty:"
    echo "${CHR_MAP}"
    exit 1
fi

if ! awk -v chr="${CHR}" '$1==chr {found=1} END{exit !found}' "${CHR_MAP}"; then
    echo "ERROR: chromosome ${CHR} is absent from chromosome map"
    exit 1
fi

########################################
# Filter and rename chromosomes
########################################

echo "========================================"
echo "Task: ${SLURM_ARRAY_TASK_ID}"
echo "Chromosome: ${CHR}"
echo "Input: ${INPUT_BCF}"
echo "Output prefix: ${OUTPUT_PREFIX}"
echo "========================================"

echo "Keeping actual biallelic SNPs and renaming chromosome..."

bcftools view \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --types snps \
    --min-alleles 2 \
    --max-alleles 2 \
    --output-type u \
    "${INPUT_BCF}" \
| bcftools annotate \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --rename-chrs "${CHR_MAP}" \
    --output-type b \
    --output "${TEMP_BCF}"

if [[ ! -s "${TEMP_BCF}" ]]; then
    echo "ERROR: temporary filtered BCF was not created"
    exit 1
fi

if ! bcftools view -h "${TEMP_BCF}" > /dev/null 2>&1; then
    echo "ERROR: temporary BCF cannot be read"
    exit 1
fi

SNP_COUNT=$(bcftools view -H "${TEMP_BCF}" | wc -l)

echo "Biallelic SNPs retained: ${SNP_COUNT}"

if [[ "${SNP_COUNT}" -eq 0 ]]; then
    echo "ERROR: no biallelic SNPs retained"
    exit 1
fi

########################################
# Convert BCF to PLINK
########################################

echo "Converting BCF to PLINK binary format..."

plink \
    --bcf "${TEMP_BCF}" \
    --double-id \
    --allow-extra-chr \
    --keep-allele-order \
    --set-missing-var-ids '@:#' \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --memory 28000 \
    --make-bed \
    --out "${OUTPUT_PREFIX}"

########################################
# Check PLINK output
########################################

for EXT in bed bim fam; do
    if [[ ! -s "${OUTPUT_PREFIX}.${EXT}" ]]; then
        echo "ERROR: PLINK output missing or empty:"
        echo "${OUTPUT_PREFIX}.${EXT}"
        exit 1
    fi
done

SAMPLE_COUNT=$(wc -l < "${OUTPUT_PREFIX}.fam")
VARIANT_COUNT=$(wc -l < "${OUTPUT_PREFIX}.bim")

echo "========================================"
echo "PLINK conversion completed"
echo "Chromosome: ${CHR}"
echo "Samples: ${SAMPLE_COUNT}"
echo "Variants: ${VARIANT_COUNT}"
echo "BED: ${OUTPUT_PREFIX}.bed"
echo "BIM: ${OUTPUT_PREFIX}.bim"
echo "FAM: ${OUTPUT_PREFIX}.fam"
echo "========================================"