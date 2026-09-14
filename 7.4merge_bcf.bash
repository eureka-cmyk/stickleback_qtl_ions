#!/bin/bash
#SBATCH --job-name=merge
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/merge-%A_%a.out


module load bcftools-uoneasy/1.19-GCC-13.2.0
vcfdir=/share/MacColl_shared/QTL/vcf/GCF_016920845.1_GAculeatus_UGA_version5_genomic/ploidy_aware_HWEPops_MQ10_BQ20_with_sex
outdir=/share/MacColl_shared/QTL/vcf/GCF_016920845.1_GAculeatus_UGA_version5_genomic/pca
mkdir -p ${outdir}

ls ${vcfdir}/*_sorted.bcf > ${outdir}/file_list.txt

echo "Merging files:"
cat ${outdir}/file_list.txt

bcftools concat \
  --naive \
  --threads 4 \
  --file-list ${outdir}/file_list.txt \
  -O b \
  -o ${outdir}/stickleback_merged_allchr.bcf

bcftools index -f ${outdir}/stickleback_merged_allchr.bcf

echo "Done!"