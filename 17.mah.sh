#!/bin/bash
#SBATCH --job-name=manhattan_6plots
#SBATCH --partition=defq
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/manhattan_6plots_%j.log

module purge
module load R-uoneasy/4.4.1-gfbf-2023b-rstudio

# -----------------------------
# Base directories
# -----------------------------
base_nine="/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_nine_main"
base_log="/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_no_outlier"

# Output directory for all 6 png files
outdir="/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_nine_main/six_main/mah_b"
mkdir -p "$outdir"

# R script path
rscript="/gpfs01/home/shxyt1/code/stickleback_qtl_ions/15.3mahn.R"

echo "Output directory: $outdir"
echo "Using R script:   $rscript"

# -----------------------------
# 4 no_outlier traits
# -----------------------------
for el in Ca_no_outlier Co_no_outlier Fe_no_outlier Sr_no_outlier; do
  echo "Running: $el"
  Rscript "$rscript" \
    --base "$base_nine" \
    --element "$el" \
    --outdir "$outdir"
done

# -----------------------------
# 2 log traits
# -----------------------------
for el in P_log V_log; do
  echo "Running: $el"
  Rscript "$rscript" \
    --base "$base_log" \
    --element "$el" \
    --outdir "$outdir"
done

echo "All 6 Manhattan plots finished."