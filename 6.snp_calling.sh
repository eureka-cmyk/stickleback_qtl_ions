#!/bin/bash
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=80g
#SBATCH --time=24:00:00
#SBATCH --job-name=call_stable
#SBATCH --array=1-22
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/call-%A_%a.out

module load bcftools-uoneasy/1.19-GCC-13.2.0
module load samtools-uoneasy/1.18-GCC-12.3.0

genome_name=GCF_016920845.1_GAculeatus_UGA_version5_genomic
reference_genome=/share/MacColl_shared/QTL/genomes/${genome_name}.fna
ploidy_file=/share/MacColl_shared/QTL/genomes/${genome_name}_ploidy.txt
GSEX_PED=/share/MacColl_shared/QTL/bams/${genome_name}/Sex_coverage/Gsex_clean.ped

chr=$(awk "NR==$SLURM_ARRAY_TASK_ID" /share/MacColl_shared/QTL/genomes/${genome_name}.bed | awk '{ print $1 }')

bamlist=/share/MacColl_shared/QTL/vcf/${genome_name}/ploidy_aware_HWEPops_MQ10_BQ20/bamlist.txt

# New output directory
outdir=/gpfs01/home/shxyt1/data/sticklebacks/BCF
mkdir -p $outdir

outprefix=${outdir}/stickleback_${chr}

echo "Calling chromosome $chr → ${outprefix}.bcf"


###############################################
# STEP 1 — mpileup (write to file, avoid pipe)
###############################################
bcftools mpileup \
  --threads 8 \
  --output-type u \
  --max-depth 5000 \
  --min-MQ 10 \
  --min-BQ 20 \
  --annotate FORMAT/DP,FORMAT/AD \
  --fasta-ref ${reference_genome} \
  --regions ${chr} \
  --bam-list ${bamlist} \
  > ${outprefix}.mpileup.bcf


###############################################
# STEP 2 — call (separate step, stable output)
###############################################
bcftools call \
  --threads 8 \
  -m \
  -P 1e-6 \
  -a GQ \
  -O b \
  --ploidy M,F \
  --ploidy-file ${ploidy_file} \
  --samples-file ${GSEX_PED} \
  -o ${outprefix}.raw.bcf \
  ${outprefix}.mpileup.bcf

rm ${outprefix}.mpileup.bcf


###############################################
# STEP 3 — sort and index
###############################################
bcftools sort \
  --temp-dir /gpfs01/home/shxyt1/temp \
  --max-mem 70G \
  -O b \
  -o ${outprefix}_sorted.bcf \
  ${outprefix}.raw.bcf

bcftools index ${outprefix}_sorted.bcf
rm ${outprefix}.raw.bcf


###############################################
# STEP 4 — Output useful QC information
###############################################
echo "--------------------------------------------"
echo "Samples in VCF:"
bcftools query -l ${outprefix}_sorted.bcf

echo "--------------------------------------------"
echo "Total variant sites:"
bcftools view -H ${outprefix}_sorted.bcf | wc -l

echo "--------------------------------------------"
echo "Total SNPs:"
bcftools view -H -v snps ${outprefix}_sorted.bcf | wc -l

echo "--------------------------------------------"
echo "Total indels:"
bcftools view -H -v indels ${outprefix}_sorted.bcf | wc -l

echo "--------------------------------------------"
echo "DONE chromosome $chr"
echo "Output: ${outprefix}_sorted.bcf"
