#!/bin/bash
#SBATCH --job-name=remontage_only
#SBATCH --partition=defq
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/remontage_only_%j.log

set -euo pipefail

module purge
module load R-uoneasy/4.4.1-gfbf-2023b-rstudio
module load ImageMagick/7.1.0-37-GCCcore-11.3.0

base="/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_log_LOD_5pc_sex"
labeled="$base/labeled_png_all"
tmp="$base/tmp_pages_from_labeled_all"
outpdf="$base/Manhattan_all_logLOD_2x3_from_labeled_png_all.pdf"

mkdir -p "$tmp"

# collect labeled pngs and sort
mapfile -t files < <(ls "$labeled"/*_labeled.png 2>/dev/null | sort)
echo "Total labeled PNGs found: ${#files[@]}"
if [[ ${#files[@]} -eq 0 ]]; then
  echo "[ABORT] No *_labeled.png found in $labeled"
  exit 1
fi

# normalize size to first image
W=$(identify -format "%w" "${files[0]}")
H=$(identify -format "%h" "${files[0]}")
echo "Target panel size: ${W}x${H}"

norm_dir="$tmp/normalized"
rm -rf "$norm_dir"
mkdir -p "$norm_dir"

norm_files=()
for f in "${files[@]}"; do
  bn=$(basename "$f")
  nf="$norm_dir/$bn"
  convert "$f" -background white -gravity center -extent "${W}x${H}" "$nf"
  norm_files+=("$nf")
done

blank="$tmp/blank.png"
convert -size "${W}x${H}" xc:white "$blank"

# clean old pages
rm -f "$tmp"/page_*.png

per_page=6
tile="2x3"
gap="+40+40"

page=1
i=0
total=${#norm_files[@]}

while [[ $i -lt $total ]]; do
  chunk=("${norm_files[@]:$i:$per_page}")
  while [[ ${#chunk[@]} -lt $per_page ]]; do
    chunk+=("$blank")
  done

  outpage="$tmp/page_${page}.png"
  echo "Creating page $page -> $outpage"

  montage "${chunk[@]}" \
    -font DejaVu-Sans \
    -tile "$tile" \
    -geometry "$gap" \
    "$outpage"

  ((page++))
  ((i+=per_page))
done

echo "Merging pages into PDF: $outpdf"
convert "$tmp"/page_*.png "$outpdf"

echo "[OK] PDF created: $outpdf"
