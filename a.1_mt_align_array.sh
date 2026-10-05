#!/bin/bash
#SBATCH --job-name=mt_align
#SBATCH --partition=defq
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/mt_align_%A_%a.log


set -euo pipefail

module load bwa-uoneasy/0.7.17-GCCcore-12.3.0
module load samtools-uoneasy/1.18-GCC-12.3.0

ref="/share/MacColl_shared/QTL/genomes/GCF_016920845.1_GAculeatus_UGA_version5_genomic.fna"
fqdir="/share/MacColl_shared/QTL/seq/Laura_Dean/HPC_data_backup/qtl/trimmed_fqs"
outdir="/gpfs01/home/shxyt1/data/sticklebacks/mt_analysis"
samples="${outdir}/samples.txt"

sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$samples")
if [[ -z "$sample" ]]; then
    echo "No sample found for array task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

tmpdir="${outdir}/tmp/${sample}"
mkdir -p "$tmpdir"

fullbam="${tmpdir}/${sample}.full.sorted.bam"
mtbam="${outdir}/${sample}.mt.bam"

trap 'rm -f "$fullbam" "${fullbam}.bai"; rm -rf "$tmpdir"' EXIT

bwa mem -t 6 \
  -R "@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA" \
  "$ref" \
  "${fqdir}/${sample}_R1.fastq.gz" \
  "${fqdir}/${sample}_R2.fastq.gz" \
  | samtools sort -@ 2 -m 4G -o "$fullbam" -

samtools index -@ 2 "$fullbam"

samtools view -@ 2 -bh -q 30 -F 0x904 \
  -o "$mtbam" "$fullbam" NC_041244.1

samtools index -@ 2 "$mtbam"

echo "Mitochondrial coverage for ${sample}:"
samtools coverage -r NC_041244.1 "$mtbam"