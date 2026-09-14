#!/bin/bash
#SBATCH --job-name=VCF_stats
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=30g
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/vcfstats-%A.out

module load bcftools-uoneasy/1.19-GCC-13.2.0

# --------------------------------------
# INPUT / OUTPUT PATHS
# --------------------------------------
VCF_DIR="/gpfs01/home/shxyt1/data/sticklebacks/BCF_clean"
OUT_DIR="/gpfs01/home/shxyt1/data/sticklebacks/BCF_clean_QC"
mkdir -p ${OUT_DIR}

echo "Running bcftools stats for all chromosomes..."
echo "VCF directory: ${VCF_DIR}"
echo "Output directory: ${OUT_DIR}"
echo "----------------------------------------------"

# --------------------------------------
# Loop through all 22 chromosomes
# --------------------------------------
for chr in {212..233}; do

    VCF="${VCF_DIR}/stickleback_NC_053${chr}.1_sorted.vcf.gz"
    STATS_OUT="${OUT_DIR}/stickleback_NC_053${chr}.1.stats"
    PLOT_DIR="${OUT_DIR}/plots_NC_053${chr}.1"

    if [[ ! -f "$VCF" ]]; then
        echo "⚠ WARNING: File not found: $VCF"
        continue
    fi
    
    echo "Processing: $VCF"
    
    # Create stats file
    bcftools stats -s - "$VCF" > "$STATS_OUT"

    # Plot results
    mkdir -p "$PLOT_DIR"
    plot-vcfstats -p "$PLOT_DIR" "$STATS_OUT"

    echo "Finished chromosome NC_053${chr}.1"
    echo "Stats: $STATS_OUT"
    echo "Plots: $PLOT_DIR/"
    echo "----------------------------------------------"

done

echo "All bcftools stats finished!"
