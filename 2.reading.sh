#!/bin/bash
# Yiwen
# 4/11/25
# for running on the UoN HPC Ada

# use devq to test code use defq for running jobs
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --array=1-383
#SBATCH --mem=64g
#SBATCH --time=48:00:00
#SBATCH --job-name=bwa_reading
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/reading-%x-%j.out

# load samtools
module load samtools-uoneasy/1.18-GCC-12.3.0
## qualimap requires java and R to be loaded
module load java-uoneasy/17.0.6
module load R-uoneasy/4.3.3-gfbf-2023b-rstudio 

# genome
genome_name=(GCF_016920845.1_GAculeatus_UGA_version5_genomic)
genome=(/share/MacColl_shared/QTL/bams/$genome_name)

#output
dir_output=(/gpfs01/home/shxyt1/data/sticklebacks/bams/stickleback_v5.0.1_assembly/)

pairdata="/gpfs01/home/shxyt1/data/sticklebacks/pheno/seq.txt"


# Debug SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_TASK_ID is: $SLURM_ARRAY_TASK_ID"

# Extract individual using awk
#individual=$(awk -F ',' 'BEGIN { OFS="," } { gsub(/^ *| *$/, "", $1); if (FNR == ENVIRON["SLURM_ARRAY_TASK_ID"]) print $1 }' $pairdata)

individual=$(awk '{print $1}' ${pairdata} | sed -n ${SLURM_ARRAY_TASK_ID}p)

# Check the result
echo "Individual extracted: $individual"

# set variables
in_filepath=($genome/raw_bams)
out_filepath=($dir_output/bamstats/QC/raw_bams)
mkdir -p $out_filepath

echo "This is job $SLURM_ARRAY_TASK_ID and will use sample $individual using bam read $in_filepath/${individual}_raw.bam"

## Create bed file to limit analysis to largest contigs (makes viewing images a pain)
## bioawk -c fastx ' length($seq) > 1000000 {print $name, 1 ,length($seq) }' $genome.bed > $genome.bed

## Run qualimap
/gpfs01/home/shxyt1/app/qualimap_v2.3/qualimap bamqc -bam $in_filepath/${individual}_raw.bam \
    -nt $SLURM_ARRAY_TASK_ID --java-mem-size=61G -outdir $out_filepath/${individual}/ -outformat HTML -outfile "${individual}_raw"
