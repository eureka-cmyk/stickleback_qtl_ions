#!/bin/bash
#SBATCH --job-name=check
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/check-%A_%a.out


module load bcftools-uoneasy/1.19-GCC-13.2.0

bcf_dir="/share/MacColl_shared/QTL/vcf/GCF_016920845.1_GAculeatus_UGA_version5_genomic/ploidy_aware_HWEPops_MQ10_BQ20_with_sex"

for bcf in ${bcf_dir}/stickleback_NC_*.bcf; do
    echo "Checking $bcf ..."
    bcftools view -H $bcf | awk '$5=="."' | wc -l
done

for bcf in ${bcf_dir}/stickleback_NC_*.bcf; do
    echo "Checking REF/ALT errors in $bcf ..."
    bcftools view -H $bcf | awk 'length($4)==0 || length($5)==0' | wc -l
done

