#!/bin/bash
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=40g
#SBATCH --time=23:00:00
#SBATCH --job-name=Stickle_call_sex_corrected
#SBATCH --array=1-22
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/snp-calling-sex-corrected-%A_%a.out

set -euo pipefail


########################################
# 1. Environment
########################################

module load bcftools-uoneasy/1.19-GCC-13.2.0
module load samtools-uoneasy/1.18-GCC-12.3.0

echo "BCFtools version:"
bcftools --version | head -1

echo "Samtools version:"
samtools --version | head -1


########################################
# 2. Reference genome
########################################

genome_name=GCF_016920845.1_GAculeatus_UGA_version5_genomic

reference_genome=/share/MacColl_shared/QTL/genomes/${genome_name}.fna

chromosome_bed=/share/MacColl_shared/QTL/genomes/${genome_name}.bed


########################################
# 3. Correct ploidy and sex files
########################################

ploidy_file=/gpfs01/home/shxyt1/data/sticklebacks/VCF/ploidy/GCF_016920845.1_v5_sex_aware_ploidy.txt

GSEX_FILE=/gpfs01/home/shxyt1/data/sticklebacks/VCF/ploidy/Gsex_bcftools.tsv


########################################
# 4. BAM list
########################################

BAMLIST=/share/MacColl_shared/QTL/vcf/${genome_name}/ploidy_aware_HWEPops_MQ10_BQ20/bamlist.txt


########################################
# 5. Quality thresholds
########################################

MQthres=10
BQthres=20


########################################
# 6. Output directory
########################################

master_output=/gpfs01/home/shxyt1/data/sticklebacks/VCF/ploidy_aware_HWEPops_MQ${MQthres}_BQ${BQthres}_sex_corrected

mkdir -p "${master_output}"


########################################
# 7. Temporary directory
########################################

TEMP_BASE=/gpfs01/home/shxyt1/temp/snp_calling_sex_corrected

TEMP_DIR=${TEMP_BASE}/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}

mkdir -p "${TEMP_DIR}"

trap 'rm -rf "${TEMP_DIR}"' EXIT


########################################
# 8. Check required input files
########################################

for input_file in \
    "${reference_genome}" \
    "${reference_genome}.fai" \
    "${chromosome_bed}" \
    "${ploidy_file}" \
    "${GSEX_FILE}" \
    "${BAMLIST}"
do
    if [[ ! -s "${input_file}" ]]; then
        echo "ERROR: required file is missing or empty:"
        echo "${input_file}"
        exit 1
    fi
done


########################################
# 9. Validate ploidy file
########################################

awk -F'\t' '
NF != 5 {
    print "ERROR: ploidy file line " NR \
          " has " NF " columns instead of 5" > "/dev/stderr"
    bad=1
}
END {
    if (NR == 0) {
        print "ERROR: ploidy file is empty" > "/dev/stderr"
        bad=1
    }
    exit bad
}
' "${ploidy_file}"

echo "Ploidy file validation: PASS"


########################################
# 10. Validate sample-sex file
########################################

awk -F'\t' '
NF != 2 {
    print "ERROR: sex file line " NR \
          " has " NF " columns instead of 2" > "/dev/stderr"
    bad=1
}

$2 != "M" && $2 != "F" {
    print "ERROR: invalid sex code at line " NR \
          ": " $2 > "/dev/stderr"
    bad=1
}

{
    sample_count[$1]++
    if (sample_count[$1] > 1) {
        print "ERROR: duplicated sample ID: " $1 > "/dev/stderr"
        bad=1
    }
}

END {
    if (NR == 0) {
        print "ERROR: sex file is empty" > "/dev/stderr"
        bad=1
    }
    exit bad
}
' "${GSEX_FILE}"

echo "Sample-sex file validation: PASS"

number_males=$(awk -F'\t' '$2=="M"{n++} END{print n+0}' "${GSEX_FILE}")

number_females=$(awk -F'\t' '$2=="F"{n++} END{print n+0}' "${GSEX_FILE}")

number_sex_samples=$(wc -l < "${GSEX_FILE}")

number_bams=$(grep -cv '^[[:space:]]*$' "${BAMLIST}")

echo "Male samples: ${number_males}"
echo "Female samples: ${number_females}"
echo "Samples in sex file: ${number_sex_samples}"
echo "BAM files in BAM list: ${number_bams}"

if [[ "${number_sex_samples}" -ne "${number_bams}" ]]; then
    echo "ERROR: number of samples in sex file does not match BAM list"
    exit 1
fi


########################################
# 11. Select chromosome
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
# 12. Output filenames
########################################

VCF_PREFIX=${master_output}/stickleback_${chr}

RAW_BCF=${VCF_PREFIX}.bcf

SORTED_BCF=${VCF_PREFIX}_sorted.bcf


########################################
# 13. Prevent accidental overwriting
########################################

if [[ -e "${SORTED_BCF}" || -e "${SORTED_BCF}.csi" ]]; then
    echo "ERROR: final output already exists:"
    echo "${SORTED_BCF}"
    echo "Delete or move it before rerunning this chromosome."
    exit 1
fi


########################################
# 14. SNP and genotype calling
########################################

echo "Starting bcftools mpileup and call..."

bcftools mpileup \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --output-type u \
    --max-depth 10000 \
    --min-MQ "${MQthres}" \
    --min-BQ "${BQthres}" \
    --platforms ILLUMINA \
    --annotate FORMAT/DP,FORMAT/AD \
    --fasta-ref "${reference_genome}" \
    --regions "${chr}" \
    --bam-list "${BAMLIST}" \
| bcftools call \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --multiallelic-caller \
    --prior 1e-6 \
    --annotate GQ \
    --output-type b \
    --ploidy-file "${ploidy_file}" \
    --samples-file "${GSEX_FILE}" \
    --output "${RAW_BCF}"

echo "Raw BCF generated:"
echo "${RAW_BCF}"

bcftools quickcheck -v "${RAW_BCF}"


########################################
# 15. Sort BCF
########################################

echo "Sorting BCF..."

bcftools sort \
    --temp-dir "${TEMP_DIR}" \
    --max-mem 25G \
    --output-type b \
    --output "${SORTED_BCF}" \
    "${RAW_BCF}"

bcftools quickcheck -v "${SORTED_BCF}"


########################################
# 16. Index sorted BCF
########################################

echo "Indexing sorted BCF..."

bcftools index \
    --force \
    --threads "${SLURM_CPUS_PER_TASK}" \
    "${SORTED_BCF}"

if [[ ! -s "${SORTED_BCF}.csi" ]]; then
    echo "ERROR: BCF index was not created"
    exit 1
fi


########################################
# 17. Remove raw unsorted BCF
########################################

rm -f "${RAW_BCF}"

echo "Removed raw unsorted BCF."


########################################
# 18. Summary
########################################

sample_total=$(bcftools query -l "${SORTED_BCF}" | wc -l)

site_total=$(bcftools index -n "${SORTED_BCF}")

snp_total=$(bcftools view -H -v snps "${SORTED_BCF}" | wc -l)

indel_total=$(bcftools view -H -v indels "${SORTED_BCF}" | wc -l)

echo "========================================"
echo "SNP calling completed successfully"
echo "Chromosome: ${chr}"
echo "Final BCF: ${SORTED_BCF}"
echo "Number of samples: ${sample_total}"
echo "Total number of sites: ${site_total}"
echo "Total number of SNPs: ${snp_total}"
echo "Total number of indels: ${indel_total}"
echo "========================================"