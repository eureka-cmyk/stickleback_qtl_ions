#!/bin/bash
#SBATCH --job-name=vcf_qc_sex_corrected
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=10g
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/qc-sex-corrected-%A_%a.out
#SBATCH --array=1-22

set -euo pipefail

module load vcftools-uoneasy/0.1.16-GCC-12.3.0

########################################
# Input and output directories
########################################

vcf_dir=/gpfs01/home/shxyt1/data/sticklebacks/VCF/ploidy_aware_HWEPops_MQ10_BQ20_sex_corrected

stat_dir=/gpfs01/home/shxyt1/data/sticklebacks/VCF/stats_sex_corrected

mkdir -p "${stat_dir}"

########################################
# Find sorted BCF files
########################################

shopt -s nullglob

bcf_files=("${vcf_dir}"/*_sorted.bcf)

if [[ ${#bcf_files[@]} -eq 0 ]]; then
    echo "ERROR: No *_sorted.bcf files found in:"
    echo "${vcf_dir}"
    exit 1
fi

echo "Number of BCF files found: ${#bcf_files[@]}"

########################################
# Select the BCF for this array task
########################################

array_index=$((SLURM_ARRAY_TASK_ID - 1))

if [[ ${array_index} -ge ${#bcf_files[@]} ]]; then
    echo "ERROR: No BCF corresponds to array task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

bcf_file="${bcf_files[$array_index]}"

if [[ ! -s "${bcf_file}" ]]; then
    echo "ERROR: BCF file is missing or empty:"
    echo "${bcf_file}"
    exit 1
fi

echo "Array task: ${SLURM_ARRAY_TASK_ID}"
echo "Processing: ${bcf_file}"

########################################
# Output prefix
########################################

base=$(basename "${bcf_file}" .bcf)

stat_prefix="${stat_dir}/${base}"

########################################
# Run VCFtools QC
########################################

echo "Running depth..."

vcftools \
    --bcf "${bcf_file}" \
    --depth \
    --out "${stat_prefix}" &

pid_depth=$!


echo "Running missing-indv..."

vcftools \
    --bcf "${bcf_file}" \
    --missing-indv \
    --out "${stat_prefix}" &

pid_missing=$!


echo "Running heterozygosity..."

vcftools \
    --bcf "${bcf_file}" \
    --het \
    --out "${stat_prefix}" &

pid_het=$!


echo "Running relatedness..."

vcftools \
    --bcf "${bcf_file}" \
    --relatedness \
    --out "${stat_prefix}" &

pid_relatedness=$!


########################################
# Wait for all analyses
########################################

wait "${pid_depth}"
wait "${pid_missing}"
wait "${pid_het}"
wait "${pid_relatedness}"

echo "========================================"
echo "QC completed successfully"
echo "Input BCF: ${bcf_file}"
echo "Results saved in: ${stat_dir}"
echo "========================================"