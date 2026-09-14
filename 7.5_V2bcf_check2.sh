#!/bin/bash
#SBATCH --job-name=check_ref_alt_V2
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=10g
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/check-ref-alt-V2-%j.out

set -euo pipefail

module load bcftools-uoneasy/1.19-GCC-13.2.0

########################################
# 1. Input and output
########################################

bcf_dir=/gpfs01/home/shxyt1/data/sticklebacks/VCF/ploidy_aware_HWEPops_MQ10_BQ20_sex_corrected

summary_file=${bcf_dir}/BCF_REF_ALT_check_V2.txt

shopt -s nullglob

bcf_files=("${bcf_dir}"/stickleback_NC_*_sorted.bcf)

if [[ ${#bcf_files[@]} -eq 0 ]]; then
    echo "ERROR: No sorted BCF files found in:"
    echo "${bcf_dir}"
    exit 1
fi

echo "Number of sorted BCF files found: ${#bcf_files[@]}"

if [[ ${#bcf_files[@]} -ne 22 ]]; then
    echo "WARNING: Expected 22 BCF files, but found ${#bcf_files[@]}"
fi

########################################
# 2. Start summary file
########################################

{
    echo -e "FILE\tSAMPLES\tTOTAL_RECORDS\tALT_DOT\tEMPTY_REF_ALT\tUNREADABLE"
} > "${summary_file}"

########################################
# 3. Check every chromosome BCF
########################################

for bcf in "${bcf_files[@]}"; do

    base=$(basename "${bcf}")

    echo "========================================"
    echo "Checking ${base}"
    echo "========================================"

    # Check whether BCF header can be read
    if ! bcftools view -h "${bcf}" > /dev/null 2>&1; then
        echo "ERROR: BCF cannot be read: ${bcf}"

        echo -e "${base}\tNA\tNA\tNA\tNA\t1" \
            >> "${summary_file}"

        continue
    fi

    samples=$(bcftools query -l "${bcf}" | wc -l)

    total_records=$(
        bcftools view -H "${bcf}" |
        wc -l
    )

    # ALT="." means no alternate allele at that site
    alt_dot=$(
        bcftools view -H "${bcf}" |
        awk '$5=="." {n++} END {print n+0}'
    )

    # Check for genuinely empty REF or ALT fields
    empty_ref_alt=$(
        bcftools view -H "${bcf}" |
        awk '
        length($4)==0 || length($5)==0 {
            n++
        }
        END {
            print n+0
        }'
    )

    echo "Samples: ${samples}"
    echo "Total records: ${total_records}"
    echo "Sites with ALT='.': ${alt_dot}"
    echo "Sites with empty REF or ALT: ${empty_ref_alt}"

    echo -e \
        "${base}\t${samples}\t${total_records}\t${alt_dot}\t${empty_ref_alt}\t0" \
        >> "${summary_file}"

done

########################################
# 4. Final report
########################################

echo "========================================"
echo "REF/ALT checking completed"
echo "Summary saved to:"
echo "${summary_file}"
echo "========================================"

column -t -s $'\t' "${summary_file}"