#!/bin/bash
#SBATCH --job-name=merge_manhattan
#SBATCH --partition=defq
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/merge_manhattan_%j.log

module purge
module load R-uoneasy/4.4.1-gfbf-2023b
module load ImageMagick/7.1.0-37-GCCcore-11.3.0


# Paths 
base="/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_nine_main"

labeled="$base/labeled_png"
tmp="$base/tmp_pages"
outpdf="$base/Manhattan_log_nine_main.pdf"

mkdir -p "$labeled" "$tmp"

echo "Base directory: $base"
echo "Output PDF: $outpdf"


# Step 1: Add centered labels
echo "=== Step 1: labeling Manhattan plots ==="

find "$base" -mindepth 2 -maxdepth 2 -type f -name "Rect_Manhtn.*_no_outlier.png" | sort | while read img; do

    element=$(basename "$img" | sed 's/^Rect_Manhtn\.//; s/\.png$//')
    out="$labeled/${element}_labeled.png"

    echo "Labeling: $element"

    convert "$img" \
        -font DejaVu-Sans \
        -gravity North \
        -pointsize 64 \
        -fill black \
        -annotate +0+30 "$element" \
        "$out"
done


# Step 2: Merge (2×4 per page)
echo "=== Step 2: merging into PDF ==="

mapfile -t files < <(ls "$labeled"/*_labeled.png | sort)

page=1
count=0
group=()

for f in "${files[@]}"; do
    group+=("$f")
    ((count++))

    if [[ $count -eq 8 ]]; then
        montage "${group[@]}" \
            -tile 2x4 \
            -geometry +30+30 \
            "$tmp/page_${page}.png"

        group=()
        count=0
        ((page++))
    fi
done

# leftover
if [[ ${#group[@]} -gt 0 ]]; then
    montage "${group[@]}" \
        -tile 2x4 \
        -geometry +30+30 \
        "$tmp/page_${page}.png"
fi


# Step 3: Convert to PDF
convert "$tmp"/page_*.png "$outpdf"

echo "=============================="
echo "PDF created successfully:"
echo "$outpdf"
echo "=============================="
