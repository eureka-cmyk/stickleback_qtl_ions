#!/bin/bash
# Yiwen
# 4/11/25
# for running on the UoN HPC Ada

# use devq to test code use defq for running jobs
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --array=1-383
#SBATCH --mem=16g
#SBATCH --time=48:00:00
#SBATCH --job-name=bwa_mapping
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/clean-%x-%j.out

# genome
genome_name=(stickleback_v5.0.1_assembly)
#/gpfs01/home/shxyt1/data/sticklebacks/bams/stickleback_v5.0.1_assembly/bamstats/QC/raw_bams/Multi-Bam-QC
bam_list="/gpfs01/home/shxyt1/data/sticklebacks/bams/$genome_name/bamstats/QC/raw_bams/Multi-Bam-QC/HiQ_bam_files.txt"

individual=$(awk "FNR==$SLURM_ARRAY_TASK_ID {print \$1}" $bam_list)
bam_file=$(awk "FNR==$SLURM_ARRAY_TASK_ID {print \$2}" $bam_list)


#master_filepath=(/share/MacColl_shared/QTL/bams/GCF_016920845.1_GAculeatus_UGA_version5_genomic)

echo "This is array task ${SLURM_ARRAY_TASK_ID}, cleaning individual $individual, using ${bam_file}"
echo "Cleaned output BAM files will be written to the folder /share/MacColl_shared/$genome_name/clean_bams"

if test -f "/share/MacColl_shared/$genome_name/clean_bams/${individual}.bam.bai"; then
    echo "${individual} already completed."
    scancel "$SLURM_JOB_ID"
else
    echo "${individual} not mapped: running bwa."
fi


module load samtools-uoneasy/1.18-GCC-12.3.0
module load bcftools-uoneasy/1.18-GCC-13.2.0

samtools view \
  --threads $SLURM_CPUS_PER_TASK \
  -q 40 \
  -f 2 \
  -F 4 \
  -b "$bam_file" | \
samtools markdup -r --threads $SLURM_CPUS_PER_TASK - \
  /share/MacColl_shared/$genome_name/clean_bams/$individual.bam

# index the final BAM files
samtools index -@ $SLURM_CPUS_PER_TASK /share/MacColl_shared/$genome_name/clean_bams/$individual.bam

echo "TEST 3"

# Test is work was completed
if test -f "/share/MacColl_shared/$genome_name/clean_bams/${individual}.bam.bai"; then
    echo "${individual} completed."
    # rm -f $master_filepath/raw_bams/${individual}_raw.bam
    # rm -f $master_filepath/raw_bams/${individual}_raw.bam.bai
else
    echo "${individual} failed to clean bam file."
fi

# check the mapping
echo "after cleaning and filtering the final mapping success was:"
samtools flagstat /share/MacColl_shared/$genome_name/clean_bams/$individual.bam