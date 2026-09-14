#!/bin/bash
#SBATCH --job-name=dist_pdf
#SBATCH --partition=defq
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/dist_pdf_%j.log

module load R-uoneasy/4.4.1-gfbf-2023b
module load ImageMagick/7.1.0-37-GCCcore-11.3.0


BASE="/gpfs01/home/shxyt1/data/sticklebacks/pheno/traits"
OUTPDF="$BASE/distribution_pdfs"
mkdir -p "$OUTPDF"

# List of transforms to process
TRANSFORMS=("original" "log" "zscore" "rint")

for T in "${TRANSFORMS[@]}"; do
    
    echo "Processing $T ..."
    
    DIR="$BASE/$T"
    TMP="$OUTPDF/tmp_$T"
    mkdir -p "$TMP"

    # 找当前变换下所有 PNG 图
    mapfile -t IMGS < <(find "$DIR" -maxdepth 1 -name "*.png" | sort)

    if [[ ${#IMGS[@]} -eq 0 ]]; then
        echo "No images found for $T, skipping."
        continue
    fi

    page=1
    group=()

    for img in "${IMGS[@]}"; do
        group+=("$img")

        if [[ ${#group[@]} -eq 16 ]]; then
            montage "${group[@]}" \
                -tile 4x4 \
                -geometry +10+10 \
                "$TMP/page_$page.png"
            group=()
            ((page++))
        fi
    done

    # leftover page
    if [[ ${#group[@]} -gt 0 ]]; then
        montage "${group[@]}" \
            -tile 4x4 \
            -geometry +10+10 \
            "$TMP/page_$page.png"
    fi

    # Convert to PDF
    convert "$TMP"/page_*.png "$OUTPDF/distribution_${T}.pdf"

    echo "Created: $OUTPDF/distribution_${T}.pdf"

done

echo "All PDFs done!"
