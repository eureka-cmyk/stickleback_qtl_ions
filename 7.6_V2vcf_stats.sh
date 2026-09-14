#!/bin/bash
#SBATCH --job-name=VCF_stats_V2
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=30g
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/vcfstats-V2-%j.out

set -euo pipefail

module load bcftools-uoneasy/1.19-GCC-13.2.0

########################################
# 1. Input and output paths
########################################

VCF_DIR=/gpfs01/home/shxyt1/data/sticklebacks/VCF/ploidy_aware_HWEPops_MQ10_BQ20_sex_corrected

OUT_DIR=/gpfs01/home/shxyt1/data/sticklebacks/VCF/bcftools_stats_sex_corrected

mkdir -p "${OUT_DIR}"

echo "BCFtools version:"
bcftools --version | head -1

echo "Input directory:"
echo "${VCF_DIR}"

echo "Output directory:"
echo "${OUT_DIR}"

echo "========================================"


########################################
# 2. Find the 22 sorted BCF files
########################################

shopt -s nullglob

BCF_FILES=(
    "${VCF_DIR}"/stickleback_NC_*_sorted.bcf
)

NUMBER_FILES=${#BCF_FILES[@]}

echo "Number of sorted BCF files found: ${NUMBER_FILES}"

if [[ "${NUMBER_FILES}" -eq 0 ]]; then
    echo "ERROR: No *_sorted.bcf files found in:"
    echo "${VCF_DIR}"
    exit 1
fi

if [[ "${NUMBER_FILES}" -ne 22 ]]; then
    echo "ERROR: Expected 22 sorted BCF files, but found ${NUMBER_FILES}"
    exit 1
fi


########################################
# 3. Run bcftools stats
########################################

for BCF in "${BCF_FILES[@]}"; do

    BASE=$(basename "${BCF}" _sorted.bcf)

    STATS_OUT=${OUT_DIR}/${BASE}.stats

    PLOT_DIR=${OUT_DIR}/plots_${BASE}

    echo "========================================"
    echo "Processing:"
    echo "${BCF}"

    # Check that the BCF can be read
    if ! bcftools view -h "${BCF}" > /dev/null 2>&1; then
        echo "ERROR: Cannot read BCF:"
        echo "${BCF}"
        exit 1
    fi

    # Check index
    if [[ ! -s "${BCF}.csi" ]]; then
        echo "ERROR: Missing BCF index:"
        echo "${BCF}.csi"
        exit 1
    fi

    echo "Running bcftools stats..."

    bcftools stats \
        --threads "${SLURM_CPUS_PER_TASK}" \
        --samples - \
        "${BCF}" \
        > "${STATS_OUT}"

    if [[ ! -s "${STATS_OUT}" ]]; then
        echo "ERROR: Stats file was not created:"
        echo "${STATS_OUT}"
        exit 1
    fi

    echo "Stats file created:"
    echo "${STATS_OUT}"


    ####################################
    # 4. Plot stats
    ####################################

    mkdir -p "${PLOT_DIR}"

    if command -v plot-vcfstats > /dev/null 2>&1; then

        echo "Generating plots..."

        plot-vcfstats \
            --prefix "${PLOT_DIR}" \
            "${STATS_OUT}"

        echo "Plots saved in:"
        echo "${PLOT_DIR}"

    else

        echo "WARNING: plot-vcfstats was not found."
        echo "The stats file was still generated successfully."

    fi

done


########################################
# 5. Create summary table
########################################

SUMMARY_FILE=${OUT_DIR}/bcftools_stats_summary_sex_corrected.tsv

printf "Chromosome\tSamples\tRecords\tSNPs\tIndels\tMultiallelic_SNPs\n" \
    > "${SUMMARY_FILE}"

for STATS_OUT in "${OUT_DIR}"/stickleback_NC_*.stats; do

    BASE=$(basename "${STATS_OUT}" .stats)

    CHR=${BASE#stickleback_}

    SAMPLES=$(
        awk -F'\t' '
        $1=="SN" && $3=="number of samples:" {
            print $4
        }
        ' "${STATS_OUT}"
    )

    RECORDS=$(
        awk -F'\t' '
        $1=="SN" && $3=="number of records:" {
            print $4
        }
        ' "${STATS_OUT}"
    )

    SNPS=$(
        awk -F'\t' '
        $1=="SN" && $3=="number of SNPs:" {
            print $4
        }
        ' "${STATS_OUT}"
    )

    INDELS=$(
        awk -F'\t' '
        $1=="SN" && $3=="number of indels:" {
            print $4
        }
        ' "${STATS_OUT}"
    )

    MULTIALLELIC=$(
        awk -F'\t' '
        $1=="SN" && $3=="number of multiallelic SNP sites:" {
            print $4
        }
        ' "${STATS_OUT}"
    )

    printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
        "${CHR}" \
        "${SAMPLES:-NA}" \
        "${RECORDS:-NA}" \
        "${SNPS:-NA}" \
        "${INDELS:-NA}" \
        "${MULTIALLELIC:-NA}" \
        >> "${SUMMARY_FILE}"

done


########################################
# 6. Final report
########################################

echo "========================================"
echo "All bcftools stats completed successfully."
echo "Stats and plots saved in:"
echo "${OUT_DIR}"
echo
echo "Summary table:"
echo "${SUMMARY_FILE}"
echo "========================================"