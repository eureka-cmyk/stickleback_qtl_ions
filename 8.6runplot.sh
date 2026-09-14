#!/bin/bash
#SBATCH --job-name=missing_and_pca_plots
#SBATCH --partition=defq
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/missing_pca_plots_%j.log

module load plink-uoneasy/1.9b_6.21-x86_64
module load python-uoneasy/3.12.3-GCCcore-13.3.0

PLINK_PREFIX="/gpfs01/home/shxyt1/data/sticklebacks/PLINK/stickleback_allchr"
OUTDIR="/gpfs01/home/shxyt1/data/sticklebacks/plink_check"

cd "${OUTDIR}"
echo "Working dir: $(pwd)"

##############################################
# [1] RUN PLINK --missing (if not already done)
##############################################

# If missingness file exists, skip calculation
if [ ! -f "${OUTDIR}/all_samples_missing.imiss" ]; then
    echo "[1] Running PLINK --missing..."
    plink \
      --bfile "${PLINK_PREFIX}" \
      --missing \
      --out "${OUTDIR}/all_samples_missing"
else
    echo "[1] all_samples_missing.imiss already exists — SKIP."
fi


##############################################
# [2] RUN PCA + MISSINGNESS PLOTS (Python)
##############################################

# Old separate plotting commands kept but disabled
# python 8.5pca_plot.py
# python 8.6missing_plot.py

echo "[2] Running combined PCA + Missingness plotting script..."
python /gpfs01/home/shxyt1/code/stickleback_qtl_ions/8.5pca_plot.py


##############################################
# [DONE]
##############################################

echo "All plots generated. Check outputs:"
echo "  ${OUTDIR}/pca_plot.png"
echo "  ${OUTDIR}/missingness_plot.png"
