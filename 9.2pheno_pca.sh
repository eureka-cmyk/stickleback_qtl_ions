#!/bin/bash
#SBATCH --job-name=pca_files
#SBATCH --partition=defq
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/pca_files_%j.log




module load plink-uoneasy/1.9b_6.21-x86_64

plink \
  --bfile /gpfs01/home/shxyt1/data/sticklebacks/PLINK/stickleback_all_clean \
  --pca 20 \
  --out   /gpfs01/home/shxyt1/data/sticklebacks/PLINK/stickleback_all_clean
