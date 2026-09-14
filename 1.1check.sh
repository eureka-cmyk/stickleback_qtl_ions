#!/bin/bash
# Yiwen
# 4/11/25
# for running on the UoN HPC Ada

# use devq to test code use defq for running jobs
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=64g
#SBATCH --time=48:00:00
#SBATCH --job-name=bwa_reading
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/reading-%x-%j.out

module load samtools-uoneasy/1.18-GCC-12.3.0
module load java-uoneasy/17.0.6
module load R-uoneasy/4.3.3-gfbf-2023b-rstudio

OUTBASE="/gpfs01/home/shxyt1/data/sticklebacks/bams/stickleback_v5.0.1_assembly/bamstats/QC/raw_bams"

for S in QTL22023 QTL22064 QTL22099 QTL22144 QTL22376; do
  BAM="/share/MacColl_shared/QTL/bams/GCF_016920845.1_GAculeatus_UGA_version5_genomic/raw_bams/${S}_raw.bam"
  OUTDIR="${OUTBASE}/${S}"
  mkdir -p "$OUTDIR"

  echo "Running qualimap for $S ..."
  ~/app/qualimap_v2.3/qualimap bamqc \
    -bam "$BAM" \
    -nt 8 \
    --java-mem-size=30G \
    -outdir "$OUTDIR" \
    -outformat HTML \
    -outfile "${S}_raw"
done
