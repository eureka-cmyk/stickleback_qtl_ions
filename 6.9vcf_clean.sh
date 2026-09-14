#!/bin/bash
#SBATCH --job-name=clean_vcf
#SBATCH --partition=defq
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/clean_vcf_%j.log

module load bcftools-uoneasy/1.19-GCC-13.2.0

bcf_dir="/share/MacColl_shared/QTL/vcf/GCF_016920845.1_GAculeatus_UGA_version5_genomic/ploidy_aware_HWEPops_MQ10_BQ20_with_sex"
outdir="/gpfs01/home/shxyt1/data/sticklebacks/BCF_clean"

mkdir -p $outdir

for bcf in ${bcf_dir}/stickleback_NC_*.bcf; do
    prefix=$(basename $bcf .bcf)
    outvcf="${outdir}/${prefix}.vcf.gz"

    echo "Cleaning $bcf → $outvcf"

    bcftools view \
        --types snps \
        --min-ac 1 \
        --exclude 'ALT="." || REF="."' \
        --force-samples \
        -Oz \
        -o $outvcf \
        $bcf

    bcftools index -t $outvcf
done

echo "All VCFs cleaned successfully!"
