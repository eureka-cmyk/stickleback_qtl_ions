#!/bin/bash
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=40g
#SBATCH --time=23:00:00
#SBATCH --job-name=sort_index_V2
#SBATCH --array=1-22
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/sort-index-V2-%A_%a.out

set -euo pipefail


########################################
# 1. Load BCFtools
########################################

module load bcftools-uoneasy/1.19-GCC-13.2.0

echo "BCFtools version:"
bcftools --version | head -1


########################################
# 2. Paths
########################################

genome_name=GCF_016920845.1_GAculeatus_UGA_version5_genomic

chromosome_bed=/share/MacColl_shared/QTL/genomes/${genome_name}.bed

bcf_dir=/gpfs01/home/shxyt1/data/sticklebacks/VCF/ploidy_aware_HWEPops_MQ10_BQ20_sex_corrected

temp_base=/gpfs01/home/shxyt1/temp/sort_index_V2

mkdir -p "${bcf_dir}"
mkdir -p "${temp_base}"


########################################
# 3. Select chromosome
########################################

chr=$(awk -v task="${SLURM_ARRAY_TASK_ID}" '
NR == task {
    print $1
    exit
}
' "${chromosome_bed}")

if [[ -z "${chr}" ]]; then
    echo "ERROR: no chromosome found for array task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

echo "========================================"
echo "Array task: ${SLURM_ARRAY_TASK_ID}"
echo "Chromosome: ${chr}"
echo "========================================"


########################################
# 4. Define input and output files
########################################

raw_bcf=${bcf_dir}/stickleback_${chr}.bcf

sorted_bcf=${bcf_dir}/stickleback_${chr}_sorted.bcf

summary_file=${bcf_dir}/stickleback_${chr}_sorting_summary.txt

temp_dir=${temp_base}/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

mkdir -p "${temp_dir}"

trap 'rm -rf "${temp_dir}"' EXIT


########################################
# 5. Handle previously completed files
########################################

if [[ -s "${sorted_bcf}" && -s "${sorted_bcf}.csi" ]]; then

    if bcftools view -h "${sorted_bcf}" > /dev/null 2>&1; then
        echo "Sorted and indexed BCF already exists:"
        echo "${sorted_bcf}"
        echo "Nothing needs to be done."
        exit 0
    fi

fi


########################################
# 6. Check raw BCF
########################################

if [[ ! -s "${raw_bcf}" ]]; then
    echo "ERROR: raw BCF is missing or empty:"
    echo "${raw_bcf}"
    exit 1
fi

echo "Checking raw BCF..."

if ! bcftools view -h "${raw_bcf}" > /dev/null; then
    echo "ERROR: raw BCF cannot be read:"
    echo "${raw_bcf}"
    exit 1
fi

echo "Raw BCF validation: PASS"


########################################
# 7. Remove incomplete previous outputs
########################################

rm -f "${sorted_bcf}"
rm -f "${sorted_bcf}.csi"


########################################
# 8. Sort BCF
########################################

echo "Sorting:"
echo "${raw_bcf}"

bcftools sort \
    --temp-dir "${temp_dir}" \
    --max-mem 20G \
    --output-type b \
    --output "${sorted_bcf}" \
    "${raw_bcf}"

if [[ ! -s "${sorted_bcf}" ]]; then
    echo "ERROR: sorted BCF was not created"
    exit 1
fi


########################################
# 9. Check sorted BCF
########################################

echo "Checking sorted BCF..."

if ! bcftools view -h "${sorted_bcf}" > /dev/null; then
    echo "ERROR: sorted BCF cannot be read:"
    echo "${sorted_bcf}"
    exit 1
fi

echo "Sorted BCF validation: PASS"


########################################
# 10. Index sorted BCF
########################################

echo "Indexing sorted BCF..."

bcftools index \
    --force \
    --threads "${SLURM_CPUS_PER_TASK}" \
    "${sorted_bcf}"

if [[ ! -s "${sorted_bcf}.csi" ]]; then
    echo "ERROR: index was not created:"
    echo "${sorted_bcf}.csi"
    exit 1
fi

echo "BCF indexing: PASS"


########################################
# 11. Generate summary statistics
########################################

echo "Calculating summary statistics..."

sample_total=$(bcftools query -l "${sorted_bcf}" | wc -l)

variant_total=$(bcftools view -H "${sorted_bcf}" | wc -l)

snp_total=$(bcftools view -H -v snps "${sorted_bcf}" | wc -l)

indel_total=$(bcftools view -H -v indels "${sorted_bcf}" | wc -l)

{
    echo "Chromosome: ${chr}"
    echo "Input raw BCF: ${raw_bcf}"
    echo "Final sorted BCF: ${sorted_bcf}"
    echo "Samples: ${sample_total}"
    echo "Total variants: ${variant_total}"
    echo "SNPs: ${snp_total}"
    echo "Indels: ${indel_total}"
} > "${summary_file}"


########################################
# 12. Delete raw BCF only after success
########################################

if [[ -s "${sorted_bcf}" && -s "${sorted_bcf}.csi" ]]; then

    bcftools view -h "${sorted_bcf}" > /dev/null

    rm -f "${raw_bcf}"

    echo "Raw unsorted BCF removed:"
    echo "${raw_bcf}"

else

    echo "ERROR: final BCF or index is missing."
    echo "Raw BCF will not be deleted."
    exit 1

fi


########################################
# 13. Final report
########################################

echo "========================================"
echo "Sorting and indexing completed"
echo "Chromosome: ${chr}"
echo "Sorted BCF: ${sorted_bcf}"
echo "Index: ${sorted_bcf}.csi"
echo "Samples: ${sample_total}"
echo "Total variants: ${variant_total}"
echo "SNPs: ${snp_total}"
echo "Indels: ${indel_total}"
echo "Summary: ${summary_file}"
echo "========================================"