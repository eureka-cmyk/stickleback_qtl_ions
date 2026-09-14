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
#SBATCH --array=1-352
#SBATCH --mem=64g
#SBATCH --time=48:00:00
#SBATCH --job-name=bwa_reading
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/reading-%x-%j.out

# load samtools
module load samtools-uoneasy/1.18-GCC-12.3.0
## qualimap requires java and R to be loaded
module load java-uoneasy/17.0.6

# Debug SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_TASK_ID is: $SLURM_ARRAY_TASK_ID"

# genome
genome_name=(GCF_016920845.1_GAculeatus_UGA_version5_genomic)


#/share/MacColl_shared/QTL/bams/GCF_016920845.1_GAculeatus_UGA_version5_genomic/clean_bams
in_filepath=(/share/MacColl_shared/QTL/bams/$genome_name/clean_bams)

# Define the pairdata file
# Input bam files
bam_list="/gpfs01/home/shxyt1/data/sticklebacks/bams/stickleback_v5.0.1_assembly/bamstats/QC/raw_bams/Multi-Bam-QC/HiQ_bam_files.txt"
# extract the individual name variable from sample name files
individual=$(awk -F ',' "FNR==$SLURM_ARRAY_TASK_ID" $bam_list | awk '{ print $1 }')
bam_file=($in_filepath/$individual.bam)


# Check the result
echo "Individual extracted: $individual"
echo "Bamfile to use is: $bam_file"

# set output variable

out_filepath=(/share/MacColl_shared/QTL/bams/$genome_name/bamstats/QC/clean_bams)
mkdir -p $out_filepath

echo "This is job $SLURM_ARRAY_TASK_ID and will use sample $individual using bam read ${bam_file}"

## Run qualimap
~/app/qualimap_v2.3/qualimap bamqc -bam ${bam_file} \
    -nt $SLURM_ARRAY_TASK_ID --java-mem-size=61G -outdir $out_filepath/${individual}/ -outformat HTML -outfile "${individual}"