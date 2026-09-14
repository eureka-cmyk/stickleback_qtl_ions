#!/bin/bash
#SBATCH --job-name=GWAS_manhattan
#SBATCH --partition=defq
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --array=1-31
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/GWAS_manhattan_%A_%a.log

module purge
module load R-uoneasy/4.4.1-gfbf-2023b-rstudio

base="/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_nine_main"

# build element list from folder names 
mapfile -t elements < <(find "$base" -maxdepth 1 -type d -name "*_no_outlier" -printf "%f\n" | sort)

total=${#elements[@]}
task=${SLURM_ARRAY_TASK_ID}

if (( task < 1 || task > total )); then
  echo "ERROR: array index $task out of range 1-$total (total=$total)"
  exit 1
fi

element="${elements[$((task-1))]}"
outdir="$base/$element"

echo "Task: $task / $total"
echo "Element: $element"
echo "Outdir:  $outdir"

Rscript code/stickleback_qtl_ions/11.3.2.manhattan.R \
  --base "$base" \
  --element "$element" \
  --outdir "$outdir"
