#!/bin/bash
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=10g
#SBATCH --time=23:00:00
#SBATCH --job-name=check-stat-V2
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/check-stat-V2-%j.out

set -euo pipefail

module load bcftools-uoneasy/1.19-GCC-13.2.0

# New sex-corrected BCF directory
vcfdir=/gpfs01/home/shxyt1/data/sticklebacks/VCF/ploidy_aware_HWEPops_MQ10_BQ20_sex_corrected

# Output summary file
summary_file=${vcfdir}/BCF_variant_statistics.txt

# Prevent an unmatched wildcard from being treated as a filename
shopt -s nullglob

bcf_files=("${vcfdir}"/*_sorted.bcf)

if [[ ${#bcf_files[@]} -eq 0 ]]; then
    echo "ERROR: No *_sorted.bcf files found in:"
    echo "${vcfdir}"
    exit 1
fi

echo "BCF directory: ${vcfdir}" > "${summary_file}"
echo "Number of BCF files: ${#bcf_files[@]}" >> "${summary_file}"
echo >> "${summary_file}"

for f in "${bcf_files[@]}"; do

    echo "Checking $(basename "${f}")"

    # Check that the BCF is readable
    if ! bcftools view -h "${f}" > /dev/null 2>&1; then
        echo "ERROR: Invalid or incomplete BCF: ${f}" | tee -a "${summary_file}"
        continue
    fi

    samples=$(bcftools query -l "${f}" | wc -l)
    total=$(bcftools view -H "${f}" | wc -l)
    snps=$(bcftools view -H -v snps "${f}" | wc -l)
    indels=$(bcftools view -H -v indels "${f}" | wc -l)

    {
        echo "===== $(basename "${f}") ====="
        echo "Samples: ${samples}"
        echo "Total variants: ${total}"
        echo "SNPs: ${snps}"
        echo "Indels: ${indels}"
        echo
    } | tee -a "${summary_file}"

done

echo "Statistics completed."
echo "Summary saved to:"
echo "${summary_file}"