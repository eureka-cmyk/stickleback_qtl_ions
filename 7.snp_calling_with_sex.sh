#!/bin/bash
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=40g
#SBATCH --time=23:00:00
#SBATCH --job-name=Stickle_call_sex
#SBATCH --array=1-22
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/snp-calling-sex-%x-%j.out

#environment
module load bcftools-uoneasy/1.19-GCC-13.2.0
module load samtools-uoneasy/1.18-GCC-12.3.0

# Genome and ploidy information
genome_name=GCF_016920845.1_GAculeatus_UGA_version5_genomic
reference_genome=/share/MacColl_shared/QTL/genomes/${genome_name}.fna
ploidy_file=/share/MacColl_shared/QTL/genomes/${genome_name}_ploidy.txt

#chromosome
chr=$(awk "NR==$SLURM_ARRAY_TASK_ID" /share/MacColl_shared/QTL/genomes/${genome_name}.bed | awk '{ print $1 }')

# Quality thresholds
MQthres=10
BQthres=20

# Master paths
master_filepath=/share/MacColl_shared/QTL
bam_dir=${master_filepath}/bams/${genome_name}/clean_bams
master_output=${master_filepath}/vcf/${genome_name}/ploidy_aware_HWEPops_MQ${MQthres}_BQ${BQthres}_with_sex
mkdir -p ${master_output}

# Output VCF/BCF prefix
VCF=stickleback_${chr}

echo "This is array task ${SLURM_ARRAY_TASK_ID}, calling SNPs for chromosome ${chr}, output: ${VCF}.bcf"

# 1) BAM list
# If you already trust the old bamlist, you can copy it; otherwise generate anew.
if [ ! -f ${master_output}/bamlist.txt ]; then
    echo "Copying existing bamlist from previous run..."
    cp /share/MacColl_shared/QTL/vcf/${genome_name}/ploidy_aware_HWEPops_MQ10_BQ20/bamlist.txt \
       ${master_output}/bamlist.txt
fi

# 2) Genomic sex file (Gsex.ped)
GSEX_PED=/share/MacColl_shared/QTL/bams/${genome_name}/Sex_coverage/Gsex_clean.ped
if [ ! -f ${GSEX_PED} ]; then
    echo "ERROR: Gsex.ped not found at ${GSEX_PED}"
    exit 1
fi

############################
# SNP and genotype calling
############################

# Generate genotype likelihoods using mpileup, then call variants with sex-aware ploidy
bcftools mpileup \
    --threads ${SLURM_CPUS_PER_TASK} \
    --output-type u \
    --max-depth 10000 \
    --min-MQ ${MQthres} \
    --min-BQ ${BQthres} \
    --platforms ILLUMINA \
    --annotate FORMAT/DP,FORMAT/AD \
    --fasta-ref ${reference_genome} \
    --regions ${chr} \
    --bam-list ${master_output}/bamlist.txt \
| bcftools call \
    --threads ${SLURM_CPUS_PER_TASK} \
    -m \
    -P 1e-6 \
    -a GQ \
    -O b \
    --ploidy M,F \
    --ploidy-file ${ploidy_file} \
    --samples-file ${GSEX_PED} \
    -o ${master_output}/${VCF}.bcf

echo "Raw BCF generated: ${master_output}/${VCF}.bcf"

############################
# Sort and index
############################

#bcftools sort -O b -o ${master_output}/${VCF}_sorted.bcf ${master_output}/${VCF}.bcf
bcftools sort \
  --temp-dir /gpfs01/home/shxyt1/temp/ \
  --max-mem 40G \
  -O b \
  -o ${master_output}/${VCF}_sorted.bcf \
  ${master_output}/${VCF}.bcf

bcftools index ${master_output}/${VCF}_sorted.bcf

# Remove unsorted BCF if sorted + index successful
if [ -f ${master_output}/${VCF}_sorted.bcf.csi ]; then
    rm ${master_output}/${VCF}.bcf
fi

############################
# Quick summary
############################

echo "These are the individuals in the VCF file:"
bcftools query -l ${master_output}/${VCF}_sorted.bcf

echo "Total number of variant sites:"
bcftools view -H ${master_output}/${VCF}_sorted.bcf | wc -l

echo "Total number of SNPs:"
bcftools view -H -v snps ${master_output}/${VCF}_sorted.bcf | wc -l

echo "Total number of indels:"
bcftools view -H -v indels ${master_output}/${VCF}_sorted.bcf | wc -l

