#!/bin/bash
#SBATCH --job-name=Manhtn_test_sexOnly
#SBATCH --partition=defq
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --array=1-3
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/Manhtn_test_sexOnly_%A_%a.log

set -euo pipefail

module purge
module load R-uoneasy/4.4.1-gfbf-2023b-rstudio

# Base directory of GWAS results
base="/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_log_LOD_5pc_sex"

# Only test these 3 elements
elements=(Ca Sr Ti)

task=${SLURM_ARRAY_TASK_ID}
element="${elements[$((task-1))]}"
outdir="${base}/${element}"

echo "Task: $task"
echo "Element: $element"
echo "Outdir: $outdir"

mlma_sex="${outdir}/${element}_MLMA_sexOnly.mlma"
mlma_std="${outdir}/${element}_MLMA.mlma"

# Make sure sexOnly MLMA exists
if [[ ! -f "$mlma_sex" ]]; then
  echo "[ERROR] Missing sexOnly MLMA file: $mlma_sex"
  exit 1
fi

# Create/refresh a symlink so the R script can always read element_MLMA.mlma
ln -sf "$(basename "$mlma_sex")" "$mlma_std"
echo "[OK] Linked: $mlma_std -> $(readlink "$mlma_std")"

# Run plotting script (unchanged)
Rscript /gpfs01/home/shxyt1/code/stickleback_qtl_ions/11.3.2.manhattan.R \
  --base "$base" \
  --element "$element" \
  --outdir "$outdir"

echo "Done: $(date)"
