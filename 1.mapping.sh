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
#SBATCH --array=1-388
#SBATCH --mem=16g
#SBATCH --time=48:00:00
#SBATCH --job-name=bwa_mapping
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/slurm-%x-%j.out

# load the necessary modules
module load bwa-uoneasy/0.7.17-GCCcore-12.3.0
module load samtools-uoneasy/1.18-GCC-12.3.0

########################
# Output and input directory

input_directory="/share/MacColl_shared/QTL/seq/Laura_Dean/HPC_data_backup/qtl/trimmed_fqs"

# Draft genome to use
genome_name="stickleback_v5.0.1_assembly"
# genome_name=(GCF_016920845.1_GAculeatus_UGA_version5_genomic)
genome="/gpfs01/home/shxyt1/data/sticklebacks/ref/fresh/$genome_name.fa"

dir_output="/gpfs01/home/shxyt1/data/sticklebacks/bams/$genome_name/raw_bams/"

/gpfs01/home/shxyt1/data/sticklebacks/bams/stickleback_v5.0.1_assembly/raw_bams/
# Make output directory
mkdir -p $dir_output

# index genome only needs doing once
## bwa index $genome

# Data on all samples
bigdata="/gpfs01/home/shxyt1/data/sticklebacks/pheno/qtl_all_phenotype_data.csv"
pairdata="/gpfs01/home/shxyt1/data/sticklebacks/pheno/seq.txt"

# awk -F "," ' FNR!=1 {print $1, "/share/MacColl_shared/QTL/seq/Laura_Dean/HPC_data_backup/qtl/trimmed_fqs/"$1"_R1.fastq.gz", "/share/MacColl_shared/QTL/seq/Laura_Dean/HPC_data_backup/qtl/trimmed_fqs/"$1"_R2.fastq.gz"}' $bigdata > $pairdata

# Gets specific sample to work with on this array
individual=$(awk "FNR == $SLURM_ARRAY_TASK_ID { print \$1 }" $pairdata)


# Use awk to process the path
# individual=$(awk -F ',' "FNR==$SLURM_ARRAY_TASK_ID" $pairdata | awk -F ',' '{ print $1 }')
forward_read=$(awk "FNR == $SLURM_ARRAY_TASK_ID { print \$2 }" $pairdata)
backward_read=$(awk "FNR == $SLURM_ARRAY_TASK_ID { print \$3 }" $pairdata)

echo "This is job $SLURM_ARRAY_TASK_ID and should use sample $individual"
echo "And should use seq files ${forward_read} and ${backward_read}"
echo "Making bam file $dir_output/${individual}_raw.bam"

## Test if file has already been created
## Once created remove raw sequence files
if test -f "$dir_output/${individual}_raw.bam.bai"; then
    echo "${individual} already completed."
    echo "This is SLURM Job ID: ${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
    scancel "${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
else
    echo "${individual} not mapped: running bwa."
fi

# Extract the header line of the fastq file
file_info=$(zcat ${forward_read} | head -n 1)

# Save the pieces of information you need as variables
flowcell_ID=$(cut -d ":" -f3 <<< "$file_info")
lane_no=$(cut -d ":" -f4 <<< "$file_info")
sample_barcode=$(cut -d ":" -f10 <<< "$file_info")

# store the read group information for the ID and PU fields as variables from the individual ones you just created
PU="${flowcell_ID}.${lane_no}.${sample_barcode}"
ID="${flowcell_ID}.${lane_no}"

bwa mem \
    -t $SLURM_CPUS_PER_TASK \
    -M \
    -v 2 \
    -R "@RG\tID:"$ID"\tSM:"$individual"\tPL:ILLUMINA\tLB:"$individual"\tPU:"$PU"" \
    $genome \
    "$forward_read" \
    "$backward_read" | \
    # mark the duplicate reads then
    # Sort the SAM files
    samtools fixmate --threads $SLURM_CPUS_PER_TASK -m -O BAM - - | \
    samtools sort --threads $SLURM_CPUS_PER_TASK -o "$dir_output/${individual}_raw.bam"

# Index bam file
samtools index "$dir_output/${individual}_raw.bam"

# Generate info  - look at how well the reads mapped
echo "The raw reads for $individual mapped with the following success:"
samtools flagstat --threads $SLURM_CPUS_PER_TASK $dir_output/${individual}_raw.bam

## Once created remove raw sequence files
if test -f "$dir_output/${individual}_raw.bam.bai"; then
    echo "${individual} successful."
else
    echo echo "${individual} unsuccessful."
fi

cd ~
