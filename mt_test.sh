#!/bin/bash
#SBATCH --job-name=mt_test
#SBATCH --partition=defq
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/mt_test_%j.log

set -euo pipefail

module load bwa-uoneasy/0.7.17-GCCcore-12.3.0
module load samtools-uoneasy/1.18-GCC-12.3.0

sample="QTL22002"
ref="/share/MacColl_shared/QTL/genomes/GCF_016920845.1_GAculeatus_UGA_version5_genomic.fna"
fqdir="/share/MacColl_shared/QTL/seq/Laura_Dean/HPC_data_backup/qtl/trimmed_fqs"
outdir="/gpfs01/home/shxyt1/data/sticklebacks/mt_test"

mkdir -p "$outdir"

bam="${outdir}/${sample}.sorted.bam"

bwa mem -t 6 \
  -R "@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA" \
  "$ref" \
  "${fqdir}/${sample}_R1.fastq.gz" \
  "${fqdir}/${sample}_R2.fastq.gz" \
  | samtools sort -@ 2 -m 4G -o "$bam" -

samtools index "$bam"

echo "Mitochondrial alignment counts:"
samtools idxstats "$bam" | awk '$1=="NC_041244.1"'

echo "Mitochondrial coverage:"
samtools coverage -r NC_041244.1 "$bam"