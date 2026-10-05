#!/bin/bash
#SBATCH --job-name=merge_noPC_manhattan
#SBATCH --partition=defq
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/merge_noPC_manhattan_%j.log

module purge
module load R-uoneasy/4.4.1-gfbf-2023b
module load ImageMagick/7.1.0-37-GCCcore-11.3.0


# ============================================================
# Paths
# ============================================================

base="/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_h_w_log_no_outlier"

plotdir="$base/plots_two_thresholds"

# ----- P < 5e-8 -----
labeled_5e8="$base/labeled_5e8"
tmp_5e8="$base/tmp_pages_5e8"

outpdf_5e8="$base/GWAS_noPC_Manhattan_threshold_5e-8_2col3row.pdf"


# ----- Bonferroni -----
labeled_bonf="$base/labeled_bonferroni"
tmp_bonf="$base/tmp_pages_bonferroni"

outpdf_bonf="$base/GWAS_noPC_Manhattan_threshold_Bonferroni_2col3row.pdf"


mkdir -p \
    "$labeled_5e8" \
    "$tmp_5e8" \
    "$labeled_bonf" \
    "$tmp_bonf"


echo "======================================================"
echo "Base directory:"
echo "$base"
echo
echo "Input Manhattan directory:"
echo "$plotdir"
echo "======================================================"



# ============================================================
# PART 1
# P < 5e-8
# ============================================================

echo
echo "======================================================"
echo "P < 5e-8"
echo "======================================================"


# ------------------------------------------------------------
# Step 1: Add label
#
# IMPORTANT:
# No resizing is performed here.
# convert only places text on the existing PNG.
# Original Manhattan aspect ratio is unchanged.
# ------------------------------------------------------------

echo "=== Step 1A: labeling P < 5e-8 Manhattan plots ==="


find "$plotdir" \
    -maxdepth 1 \
    -type f \
    -name "Rect_Manhtn.*_threshold_5e-8.png" \
    | sort \
    | while read img
do

    filename=$(basename "$img")

    element=$(echo "$filename" \
        | sed 's/^Rect_Manhtn\.//; s/_threshold_5e-8\.png$//')

    out="$labeled_5e8/${element}_labeled.png"

    echo "Labeling: $element"

    convert "$img" \
        -font DejaVu-Sans \
        -gravity North \
        -pointsize 64 \
        -fill black \
        -annotate +0+30 "${element}   P < 5e-8" \
        "$out"

done



# ------------------------------------------------------------
# Step 2: Merge 6 plots per page
# 2 columns × 3 rows
#
# IMPORTANT:
# -geometry +30+30 contains NO width/height.
# Therefore montage does NOT resize the Manhattan plots.
# It only adds spacing between them.
# ------------------------------------------------------------

echo
echo "=== Step 2A: merging P < 5e-8 plots (2 x 3 per page) ==="


mapfile -t files_5e8 < <(
    find "$labeled_5e8" \
        -maxdepth 1 \
        -type f \
        -name "*_labeled.png" \
        | sort
)


echo "Number of plots: ${#files_5e8[@]}"


page=1
count=0
group=()


for f in "${files_5e8[@]}"; do

    group+=("$f")

    count=$((count + 1))


    if [[ $count -eq 6 ]]; then

        echo "Creating page $page"

        montage "${group[@]}" \
            -tile 2x3 \
            -geometry +30+30 \
            -background white \
            "$tmp_5e8/page_${page}.png"

        group=()
        count=0
        page=$((page + 1))

    fi

done


# leftover plots
if [[ ${#group[@]} -gt 0 ]]; then

    echo "Creating final page $page with ${#group[@]} plots"

    montage "${group[@]}" \
        -tile 2x3 \
        -geometry +30+30 \
        -background white \
        "$tmp_5e8/page_${page}.png"

fi



# ------------------------------------------------------------
# Step 3: Convert pages to PDF
# ------------------------------------------------------------

echo
echo "=== Step 3A: creating P < 5e-8 PDF ==="


convert \
    "$tmp_5e8"/page_*.png \
    "$outpdf_5e8"


echo
echo "P < 5e-8 PDF created:"
echo "$outpdf_5e8"



# ============================================================
# PART 2
# Bonferroni
# ============================================================

echo
echo "======================================================"
echo "Bonferroni threshold"
echo "P < 1.099e-8"
echo "======================================================"


# ------------------------------------------------------------
# Step 1: Add label
# ------------------------------------------------------------

echo "=== Step 1B: labeling Bonferroni Manhattan plots ==="


find "$plotdir" \
    -maxdepth 1 \
    -type f \
    -name "Rect_Manhtn.*_threshold_Bonferroni.png" \
    | sort \
    | while read img
do

    filename=$(basename "$img")

    element=$(echo "$filename" \
        | sed 's/^Rect_Manhtn\.//; s/_threshold_Bonferroni\.png$//')

    out="$labeled_bonf/${element}_labeled.png"

    echo "Labeling: $element"

    convert "$img" \
        -font DejaVu-Sans \
        -gravity North \
        -pointsize 64 \
        -fill black \
        -annotate +0+30 "${element}   P < 1.099e-8" \
        "$out"

done



# ------------------------------------------------------------
# Step 2: Merge 6 plots per page
# 2 columns × 3 rows
# ------------------------------------------------------------

echo
echo "=== Step 2B: merging Bonferroni plots (2 x 3 per page) ==="


mapfile -t files_bonf < <(
    find "$labeled_bonf" \
        -maxdepth 1 \
        -type f \
        -name "*_labeled.png" \
        | sort
)


echo "Number of plots: ${#files_bonf[@]}"


page=1
count=0
group=()


for f in "${files_bonf[@]}"; do

    group+=("$f")

    count=$((count + 1))


    if [[ $count -eq 6 ]]; then

        echo "Creating page $page"

        montage "${group[@]}" \
            -tile 2x3 \
            -geometry +30+30 \
            -background white \
            "$tmp_bonf/page_${page}.png"

        group=()
        count=0
        page=$((page + 1))

    fi

done


# leftover plots
if [[ ${#group[@]} -gt 0 ]]; then

    echo "Creating final page $page with ${#group[@]} plots"

    montage "${group[@]}" \
        -tile 2x3 \
        -geometry +30+30 \
        -background white \
        "$tmp_bonf/page_${page}.png"

fi



# ------------------------------------------------------------
# Step 3: Convert pages to PDF
# ------------------------------------------------------------

echo
echo "=== Step 3B: creating Bonferroni PDF ==="


convert \
    "$tmp_bonf"/page_*.png \
    "$outpdf_bonf"


echo
echo "Bonferroni PDF created:"
echo "$outpdf_bonf"



# ============================================================
# Finished
# ============================================================

echo
echo "======================================================"
echo "ALL DONE"
echo "======================================================"

echo
echo "P < 5e-8:"
echo "$outpdf_5e8"

echo
echo "Bonferroni P < 1.099e-8:"
echo "$outpdf_bonf"

echo
echo "Layout:"
echo "2 columns x 3 rows per page"

echo
echo "IMPORTANT:"
echo "Original Manhattan plot dimensions/aspect ratios were preserved."
echo "No resizing was applied during montage."

echo "======================================================"