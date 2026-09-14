#!/bin/bash
#SBATCH --job-name=merge_manhattan_sig6
#SBATCH --partition=defq
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/merge_manhattan_sig6_%j.log

module purge
module load R-uoneasy/4.4.1-gfbf-2023b
module load ImageMagick/7.1.0-37-GCCcore-11.3.0

# ------------------------------------------------------------
# Base directories
# ------------------------------------------------------------
base_nine="/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_nine_main"
base_log="/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_no_outlier"

labeled="$base_nine/labeled_png_sig6"
tmp="$base_nine/tmp_pages_sig6"
outpdf="$base_nine/Manhattan_sig6.pdf"

mkdir -p "$labeled" "$tmp"

# ------------------------------------------------------------
# Elements
# ------------------------------------------------------------
elements_no_outlier=(Ca Co Fe Sr)
elements_log=(P V)

# layout order (important)
layout=(Ca Co Fe Sr P V)

echo "Output PDF: $outpdf"

# ------------------------------------------------------------
# Step 1: Add LARGE titles (element name only)
# ------------------------------------------------------------
files=()
missing=0

for el in "${layout[@]}"; do

    if [[ " ${elements_no_outlier[*]} " =~ " ${el} " ]]; then
        img="$base_nine/${el}_no_outlier/Rect_Manhtn.${el}_no_outlier.png"
    else
        img="$base_log/${el}_log/Rect_Manhtn.${el}_log.png"
    fi

    if [[ ! -f "$img" ]]; then
        echo "[ERROR] Missing Manhattan plot for $el:"
        echo "        $img"
        ((missing++))
        continue
    fi

    out="$labeled/${el}_labeled.png"
    echo "Labeling: $el"

    convert "$img" \
        -font DejaVu-Sans-Bold \
        -gravity North \
        -pointsize 96 \
        -fill black \
        -annotate +0+40 "$el" \
        "$out"

    files+=("$out")
done

if [[ $missing -gt 0 ]]; then
    echo "[ABORT] $missing plot(s) missing. Fix paths before rerun."
    exit 1
fi

# ------------------------------------------------------------
# Step 2: Merge into 1 page (2 columns x 3 rows)
# ------------------------------------------------------------
montage "${files[@]}" \
    -tile 2x3 \
    -geometry +40+40 \
    "$tmp/page_1.png"

# ------------------------------------------------------------
# Step 3: Convert to PDF
# ------------------------------------------------------------
convert "$tmp/page_1.png" "$outpdf"

echo "=============================="
echo "PDF created successfully:"
echo "$outpdf"
echo "=============================="
