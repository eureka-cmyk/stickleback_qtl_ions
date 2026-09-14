#!/bin/bash
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=10g
#SBATCH --time=23:00:00
#SBATCH --job-name=check-stat
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/check-stat-%x-%j.out

vcfdir=/share/MacColl_shared/QTL/vcf/GCF_016920845.1_GAculeatus_UGA_version5_genomic/ploidy_aware_HWEPops_MQ10_BQ20_with_sex
module load bcftools-uoneasy/1.19-GCC-13.2.0

for f in ${vcfdir}/*_sorted.bcf; do
    echo "===== $(basename $f) ====="
    echo "Samples:"
    bcftools query -l $f | wc -l

    echo "Total variants:"
    bcftools view -H $f | wc -l

    echo "SNPs:"
    bcftools view -H -v snps $f | wc -l

    echo "Indels:"
    bcftools view -H -v indels $f | wc -l

    echo
done
