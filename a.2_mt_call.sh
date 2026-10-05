#!/bin/bash
#SBATCH --job-name=mt_call
#SBATCH --partition=defq
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/mt_call_%j.log

set -euo pipefail

module load bcftools-uoneasy/1.19-GCC-13.2.0

outdir="/gpfs01/home/shxyt1/data/sticklebacks/mt_analysis"
ref="/share/MacColl_shared/QTL/genomes/GCF_016920845.1_GAculeatus_UGA_version5_genomic.fna"

bcftools mpileup \
  -f "$ref" \
  -r NC_041244.1 \
  -b "$outdir/mt_bams.list" \
  -q 30 -Q 20 -d 10000 -D \
  -a FORMAT/DP,FORMAT/AD \
  -Ou \
  | bcftools call \
      -S "$outdir/haploid_samples.txt" \
      -m -v --threads 3 \
      -Oz -o "$outdir/mitochondrial_variants.vcf.gz"

bcftools index -t "$outdir/mitochondrial_variants.vcf.gz"