#!/bin/bash
#SBATCH --job-name=bcftovcf
#SBATCH --partition=defq
#SBATCH --cpus-per-task=4
#SBATCH --array=1-22
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/plink_convert_%j.log

#environment
module load bcftools-uoneasy/1.19-GCC-13.2.0

bcf_dir="/share/MacColl_shared/QTL/vcf/GCF_016920845.1_GAculeatus_UGA_version5_genomic/ploidy_aware_HWEPops_MQ10_BQ20_with_sex"
vcf_dir=/gpfs01/home/shxyt1/data/sticklebacks/VCF

bcf_list=($(ls ${bcf_dir}/stickleback_NC_*.bcf | sort))

this_bcf=${bcf_list[$SLURM_ARRAY_TASK_ID-1]}
prefix=$(basename $this_bcf .bcf)

echo "Converting: $this_bcf → ${prefix}.vcf.gz"

bcftools view -Ov $this_bcf | bgzip -@ 4 > ${vcf_dir}/${prefix}.vcf.gz

bcftools index -t ${vcf_dir}/${prefix}.vcf.gz

echo "Done: ${prefix}.vcf.gz"
