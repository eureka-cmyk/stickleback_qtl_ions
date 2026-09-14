#!/bin/bash
#SBATCH --job-name=fix_missing_ids
#SBATCH --partition=defq
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/fix_missing_ids_%j.log

PLINK_DIR="/gpfs01/home/shxyt1/data/sticklebacks/PLINK"
cd $PLINK_DIR

echo "Fixing missing SNP IDs (. → CHR:BP)..."

for bim in stickleback_NC_*.bim; do
    echo "Processing $bim"
    awk 'BEGIN{OFS="\t"} {
        if ($2=="." || $2=="") $2 = $1 ":" $4;
        print
    }' $bim > ${bim}.tmp
    mv ${bim}.tmp $bim
done

echo "DONE — SNP IDs fixed."
