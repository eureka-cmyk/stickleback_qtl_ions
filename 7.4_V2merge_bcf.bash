#!/bin/bash
#SBATCH --job-name=merge_sex_corrected
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=40g
#SBATCH --time=23:00:00
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/merge-sex-corrected-%j.out

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

# New sex-corrected chromosome BCF files
vcfdir=/gpfs01/home/shxyt1/data/sticklebacks/VCF/ploidy_aware_HWEPops_MQ10_BQ20_sex_corrected

# New output directory
outdir=/gpfs01/home/shxyt1/data/sticklebacks/VCF/pca_sex_corrected

mkdir -p "${outdir}"


########################################
# 3. Output files
########################################

file_list=${outdir}/file_list_sex_corrected.txt

merged_bcf=${outdir}/stickleback_merged_allchr_sex_corrected.bcf


########################################
# 4. Build file list in chromosome order
########################################

> "${file_list}"

while read -r chr rest; do

    bcf_file=${vcfdir}/stickleback_${chr}_sorted.bcf

    if [[ ! -s "${bcf_file}" ]]; then
        echo "ERROR: BCF file is missing or empty:"
        echo "${bcf_file}"
        exit 1
    fi

    if [[ ! -s "${bcf_file}.csi" ]]; then
        echo "ERROR: BCF index is missing:"
        echo "${bcf_file}.csi"
        exit 1
    fi

    echo "${bcf_file}" >> "${file_list}"

done < <(head -22 "${chromosome_bed}")


########################################
# 5. Check file list
########################################

number_files=$(wc -l < "${file_list}")

echo "Number of chromosome BCF files: ${number_files}"

if [[ "${number_files}" -ne 22 ]]; then
    echo "ERROR: expected 22 chromosome BCF files, but found ${number_files}"
    exit 1
fi

echo "Files will be merged in this order:"
cat "${file_list}"


########################################
# 6. Prevent accidental overwriting
########################################

if [[ -e "${merged_bcf}" || -e "${merged_bcf}.csi" ]]; then
    echo "ERROR: merged output already exists:"
    echo "${merged_bcf}"
    echo "Move or delete the existing file before rerunning."
    exit 1
fi


########################################
# 7. Concatenate chromosomes
########################################

echo "Starting BCF concatenation..."

bcftools concat \
    --naive \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --file-list "${file_list}" \
    --output-type b \
    --output "${merged_bcf}"


########################################
# 8. Check merged BCF
########################################

if [[ ! -s "${merged_bcf}" ]]; then
    echo "ERROR: merged BCF was not created."
    exit 1
fi

if ! bcftools view -h "${merged_bcf}" > /dev/null; then
    echo "ERROR: merged BCF cannot be read."
    exit 1
fi

echo "Merged BCF validation: PASS"


########################################
# 9. Index merged BCF
########################################

echo "Indexing merged BCF..."

bcftools index \
    --force \
    --threads "${SLURM_CPUS_PER_TASK}" \
    "${merged_bcf}"

if [[ ! -s "${merged_bcf}.csi" ]]; then
    echo "ERROR: merged BCF index was not created."
    exit 1
fi


########################################
# 10. Summary
########################################

sample_total=$(bcftools query -l "${merged_bcf}" | wc -l)

contig_total=$(bcftools query -f '%CHROM\n' "${merged_bcf}" | uniq | wc -l)

echo "========================================"
echo "BCF concatenation completed successfully"
echo "Merged BCF: ${merged_bcf}"
echo "Index: ${merged_bcf}.csi"
echo "Samples: ${sample_total}"
echo "Chromosomes represented: ${contig_total}"
echo "========================================"