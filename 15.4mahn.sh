#!/bin/bash
#SBATCH --job-name=GWAS_manhattan_logLOD
#SBATCH --partition=defq
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --array=1-31
#SBATCH --time=06:00:00
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/GWAS_manhattan_logLOD_%A_%a.log

module purge
module load R-uoneasy/4.4.1-gfbf-2023b-rstudio

# Base directory of log_LOD GWAS results
base="/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_log_LOD_5pc_sex"

# Build element list from folder names (exclude GRM)
mapfile -t elements < <(find "$base" -maxdepth 1 -mindepth 1 -type d -printf "%f\n" \
  | grep -v '^GRM$' \
  | sort)

total=${#elements[@]}
task=${SLURM_ARRAY_TASK_ID}

echo "Total element folders found: $total"

if (( task < 1 || task > total )); then
  echo "ERROR: array index $task out of range 1-$total (total=$total)"
  exit 1
fi

element="${elements[$((task-1))]}"
outdir="$base/$element"

echo "Task: $task / $total"
echo "Element: $element"
echo "Outdir:  $outdir"

# Run plotting script
Rscript /gpfs01/home/shxyt1/code/stickleback_qtl_ions/11.3.2.manhattan.R \
  --base "$base" \
  --element "$element" \
  --outdir "$outdir"
