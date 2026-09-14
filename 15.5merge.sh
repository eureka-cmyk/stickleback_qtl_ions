#!/bin/bash
#SBATCH --job-name=merge_manhattan_all_logLOD
#SBATCH --partition=defq
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/merge_manhattan_all_logLOD_%j.log

module purge
module unload ImageMagick
module load R-uoneasy/4.4.1-gfbf-2023b


# ------------------------------------------------------------
# Base directory (LOG_LOD GWAS results)
# ------------------------------------------------------------
base="/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_log_LOD_5pc_sex"

# Output folders/files
labeled="$base/labeled_png_all"
tmp="$base/tmp_pages_all"
outpdf="$base/Manhattan_all_logLOD.pdf"

# Output tables
out_table="$base/All_log_LOD_TopSNPs_sig_sugg.txt"
table_png="$base/All_log_LOD_TopSNPs_sig_sugg.png"
table_pdf="$base/All_log_LOD_TopSNPs_sig_sugg.pdf"

mkdir -p "$labeled" "$tmp"

# ------------------------------------------------------------
# Thresholds
# ------------------------------------------------------------
sig=5e-8
sugg=1e-6

# ------------------------------------------------------------
# Collect element directories (exclude GRM and helper folders)
# ------------------------------------------------------------
mapfile -t elements < <(find "$base" -maxdepth 1 -mindepth 1 -type d -printf "%f\n" \
  | grep -v '^GRM$' \
  | grep -v '^labeled_png_' \
  | grep -v '^tmp_pages_' \
  | sort)

total=${#elements[@]}

echo "Base:   $base"
echo "Total element dirs: $total"
echo "Output PDF: $outpdf"
echo "Output table: $out_table"

if [[ $total -lt 1 ]]; then
  echo "[ABORT] No element folders found under: $base"
  exit 1
fi

# ------------------------------------------------------------
# Step A: Build the significant/suggestive table from MLMA
# Output columns:
# Element  SNP  CHR  BP  P  Category
# ------------------------------------------------------------
echo -e "Element\tSNP\tCHR\tBP\tP\tCategory" > "$out_table"

for el in "${elements[@]}"; do
  mlma="$base/$el/${el}_MLMA.mlma"

  if [[ ! -f "$mlma" ]]; then
    echo -e "${el}\tNA\tNA\tNA\tNA\tmissing_MLMA" >> "$out_table"
    continue
  fi

  # Robust AWK: detect column indices from header row
  awk -v el="$el" -v sig="$sig" -v sugg="$sugg" '
    BEGIN{
      hits=0
      idxSNP=idxCHR=idxBP=idxP=0
    }
    NR==1{
      for(i=1;i<=NF;i++){
        if($i=="SNP") idxSNP=i
        if($i=="Chr" || $i=="CHR") idxCHR=i
        if($i=="bp"  || $i=="BP")  idxBP=i
        if($i=="p"   || $i=="P")   idxP=i
      }
      next
    }
    NR>1{
      if(idxSNP==0 || idxCHR==0 || idxBP==0 || idxP==0) next
      snp=$idxSNP; chr=$idxCHR; bp=$idxBP; p=$idxP

      # skip missing
      if(p=="NA" || p=="" || chr=="NA" || bp=="NA") next

      # numeric compare
      if(p+0 < sig+0){
        printf("%s\t%s\t%s\t%s\t%.10g\tsignificant\n", el, snp, chr, bp, p)
        hits=1
      } else if(p+0 < sugg+0){
        printf("%s\t%s\t%s\t%s\t%.10g\tsuggestive\n", el, snp, chr, bp, p)
        hits=1
      }
    }
    END{
      if(hits==0){
        printf("%s\tNA\tNA\tNA\tNA\tno_hit\n", el)
      }
    }
  ' "$mlma" >> "$out_table"
done

echo "[OK] Table saved: $out_table"

# ------------------------------------------------------------
# Step B: Label each Manhattan PNG with element name (big title)
# Expect the manhattan png name:  <element>.png  inside element folder
# (Produced by CMplot with file.name = element)
# ------------------------------------------------------------
files=()
missing=0

for el in "${elements[@]}"; do
  # try standard CMplot name first
  img="$base/$el/${el}.png"

# if not found, try Rect_Manhtn.* pattern
  if [[ ! -f "$img" ]]; then
      img_alt=$(ls "$base/$el"/Rect_Manhtn*.png 2>/dev/null | head -n 1)
      if [[ -n "$img_alt" ]]; then
          img="$img_alt"
      fi
  fi


  out="$labeled/${el}_labeled.png"
  echo "Labeling: $el"

  convert "$img" \
    -font DejaVu-Sans-Bold \
    -gravity North \
    -pointsize 72 \
    -fill black \
    -annotate +0+30 "$el" \
    "$out"

  files+=("$out")
done

if [[ $missing -gt 0 ]]; then
  echo "[ABORT] $missing plot(s) missing. Fix paths before rerun."
  exit 1
fi

echo "[OK] Labeled images: ${#files[@]}"

# ------------------------------------------------------------
# Step C: Montage into multiple pages (fixed grid)
# Here: 3x3 per page (9 plots) to keep sizes identical
# ------------------------------------------------------------
tile="3x3"
per_page=9
gap="+30+30"

page=1
i=0
n=${#files[@]}

while [[ $i -lt $n ]]; do
  chunk=("${files[@]:$i:$per_page}")

  # Fill empty slots to keep identical panel sizes on the last page
  if [[ ${#chunk[@]} -lt $per_page ]]; then
    # Use a blank spacer based on the first image size
    blank="$tmp/_blank.png"
    if [[ ! -f "$blank" ]]; then
      convert "${files[0]}" -alpha off -fill white -colorize 100% "$blank"
    fi
    while [[ ${#chunk[@]} -lt $per_page ]]; do
      chunk+=("$blank")
    done
  fi

  outpage="$tmp/page_${page}.png"
  echo "Creating montage page $page -> $outpage"

  montage "${chunk[@]}" \
  -font DejaVu-Sans \
  -tile "$tile" \
  -geometry "$gap" \
  "$outpage"


  ((page++))
  ((i+=per_page))
done

# ------------------------------------------------------------
# Step D: Render the table into a PNG + PDF (1 page)
# Use R to make a clean table image (top N rows only to keep readable)
# ------------------------------------------------------------
Rscript --vanilla - <<'RSCRIPT'
library(data.table)
library(gridExtra)
library(grid)

base <- "/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_log_LOD_5pc_sex"
infile <- file.path(base, "All_log_LOD_TopSNPs_sig_sugg.txt")
pngfile <- file.path(base, "All_log_LOD_TopSNPs_sig_sugg.png")

dt <- fread(infile)

# Keep only significant/suggestive rows; keep "no_hit" as well (optional)
# If you want ONLY sig+sugg, uncomment the next line:
# dt <- dt[Category %in% c("significant","suggestive")]

# If too many rows, show first 200 to keep the PDF readable
max_rows <- 200
dt_show <- dt
if (nrow(dt_show) > max_rows) dt_show <- dt_show[1:max_rows]

# Render as a table to PNG
png(pngfile, width = 2400, height = 3200, res = 200)
grid.newpage()
title <- textGrob("Significant & Suggestive SNPs (log_LOD GWAS)", gp=gpar(fontsize=16, fontface="bold"))
tbl <- tableGrob(dt_show, rows = NULL, theme = ttheme_minimal(base_size = 8))
grid.arrange(title, tbl, ncol=1, heights=c(0.06, 0.94))
dev.off()
RSCRIPT

# Convert table png to pdf (single page)
convert "$table_png" "$table_pdf"
echo "[OK] Table PDF saved: $table_pdf"

# ------------------------------------------------------------
# Step E: Convert all montage pages + table page into one PDF
# ------------------------------------------------------------
echo "Merging pages into final PDF: $outpdf"
convert "$tmp"/page_*.png "$table_png" "$outpdf"

echo "=============================="
echo "Final PDF created:"
echo "$outpdf"
echo "Table file:"
echo "$out_table"
echo "=============================="
