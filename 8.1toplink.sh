#!/bin/bash
#SBATCH --job-name=plink_vcf2bed
#SBATCH --partition=defq
#SBATCH --cpus-per-task=4
#SBATCH --array=1-22
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/plink_vcf2bed_%A_%a.log

module load bcftools-uoneasy/1.19-GCC-13.2.0
module load plink-uoneasy/1.9b_6.21-x86_64

vcf_dir="/gpfs01/home/shxyt1/data/sticklebacks/BCF_clean"
outdir="/gpfs01/home/shxyt1/data/sticklebacks/PLINK"
mkdir -p $outdir

vcf_list=($(ls ${vcf_dir}/stickleback_NC_*.vcf.gz | sort))
this_vcf=${vcf_list[$SLURM_ARRAY_TASK_ID-1]}
prefix=$(basename $this_vcf .vcf.gz)

echo "Processing $this_vcf → PLINK (with chromosome rename)"

bcftools annotate \
    --rename-chrs /gpfs01/home/shxyt1/data/sticklebacks/chr_map.txt \
    $this_vcf \
| plink \
    --vcf /dev/stdin \
    --double-id \
    --allow-extra-chr \
    --make-bed \
    --out ${outdir}/${prefix}

echo "Done: ${prefix}.bed"
