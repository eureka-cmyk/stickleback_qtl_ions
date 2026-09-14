#!/bin/bash
#SBATCH --job-name=plink_clean
#SBATCH --partition=defq
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/plink_clean_%j.log

module load plink-uoneasy/1.9b_6.21-x86_64


plink \
  --bfile /gpfs01/home/shxyt1/data/sticklebacks/PLINK/stickleback_allchr \
  --geno 0.05 \
  --maf 0.01 \
  --make-bed \
  --out /gpfs01/home/shxyt1/data/sticklebacks/PLINK/stickleback_all_clean
