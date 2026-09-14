#!/bin/bash
#SBATCH --job-name=manhattan_6plots_EDA10
#SBATCH --partition=defq
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/manhattan_6plots_EDA10_%j.log

module purge
module load R-uoneasy/4.4.1-gfbf-2023b-rstudio

# -----------------------------
# Base directory
# -----------------------------
base="/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_EDA10"

# Output directory for all 6 png files
outdir="/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_EDA10/mah_b"
mkdir -p "$outdir"

# R script path
rscript="/gpfs01/home/shxyt1/code/stickleback_qtl_ions/15.3mahn.R"

echo "Base directory:   $base"
echo "Output directory: $outdir"
echo "Using R script:   $rscript"

# -----------------------------
# Six main traits
# -----------------------------
elements=(
  Ca_no_outlier
  Co_no_outlier
  Fe_no_outlier
  Sr_no_outlier
  P_log
  V_log
)

for el in "${elements[@]}"; do
  echo "=============================="
  echo "Running Manhattan plot for: $el"
  echo "Input file:"
  echo "$base/$el/${el}_MLMA.mlma"
  echo "=============================="

  Rscript "$rscript" \
    --base "$base" \
    --element "$el" \
    --outdir "$outdir"
done

echo "=============================="
echo "All 6 EDA10 Manhattan plots finished."
echo "Output directory:"
echo "$outdir"
echo "=============================="