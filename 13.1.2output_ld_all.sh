#!/bin/bash
#SBATCH --job-name=ld_lead_fullwindow
#SBATCH --partition=defq
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/ld_lead_fullwindow_%j.log
module load plink-uoneasy/1.9b_6.21-x86_64

bfile=/gpfs01/home/shxyt1/data/sticklebacks/PLINK/stickleback_all_clean
in=/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_EDA10/locus_tbl2_for_ld.tsv
#in=/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_nine_main/six_main/locus_tbl2_for_ld.tsv

outdir=/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_nine_main/six_main

mkdir -p "$outdir"

tail -n +2 "$in" | while IFS=$'\t' read -r Element CHR locus_id lead_SNP window_start window_end; do
  outprefix="${outdir}/${Element}_chr${CHR}_locus${locus_id}_lead${lead_SNP}"

  # 1) snplist for all SNPs in the window (from the bfile)
  plink --bfile "$bfile" \
    --chr "$CHR" \
    --from-bp "$window_start" \
    --to-bp "$window_end" \
    --write-snplist \
    --out "$outprefix"

  # 2) LD: lead SNP vs ALL window SNPs (no truncation)
  plink --bfile "$bfile" \
    --extract "${outprefix}.snplist" \
    --r2 \
    --ld-snp "$lead_SNP" \
    --ld-window 999999 \
    --ld-window-kb 1000000 \
    --ld-window-r2 0 \
    --out "$outprefix"

done



#!/bin/bash
#SBATCH --job-name=LD_from_locus_table
#SBATCH --partition=defq
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/LD_from_locus_table_%j.log

set -euo pipefail

echo "Start time: $(date)"
echo "Running on: $(hostname)"

# ==============================
# 1. Paths
# ==============================
PLINK_BIN="plink"
BFILE="/gpfs01/home/shxyt1/data/sticklebacks/PLINK/stickleback_all_clean"
LOCUS_TBL="/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_nine_main/six_main/locus_tbl2_for_ld.tsv"
OUTDIR="/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_nine_main/six_main"

mkdir -p "${OUTDIR}"

# ==============================
# 2. Check input files
# ==============================
if [[ ! -f "${BFILE}.bed" ]]; then
    echo "Error: ${BFILE}.bed not found"
    exit 1
fi

if [[ ! -f "${BFILE}.bim" ]]; then
    echo "Error: ${BFILE}.bim not found"
    exit 1
fi

if [[ ! -f "${BFILE}.fam" ]]; then
    echo "Error: ${BFILE}.fam not found"
    exit 1
fi

if [[ ! -f "${LOCUS_TBL}" ]]; then
    echo "Error: ${LOCUS_TBL} not found"
    exit 1
fi

# ==============================
# 3. Run LD for each locus
# ==============================
tail -n +2 "${LOCUS_TBL}" | \
while IFS=$'\t' read -r Element CHR locus_id n_snps minBP maxBP spanBP lead_SNP lead_BP lead_P flank_kb window_start window_end
do
    prefix="${OUTDIR}/${Element}_chr${CHR}_locus${locus_id}"

    echo "----------------------------------------"
    echo "Processing: Element=${Element}, CHR=${CHR}, locus_id=${locus_id}"
    echo "Lead SNP: ${lead_SNP}"
    echo "Window: ${window_start}-${window_end}"

    # Generate SNP list from BIM file based on chromosome and BP range
    awk -v chr="${CHR}" -v start="${window_start}" -v end="${window_end}" '
        BEGIN {FS=OFS="\t"}
        $1 == chr && $4 >= start && $4 <= end {print $2}
    ' "${BFILE}.bim" > "${prefix}.snplist"

    # Check whether snplist is empty
    if [[ ! -s "${prefix}.snplist" ]]; then
        echo "Warning: ${prefix}.snplist is empty, skipping this locus."
        continue
    fi

    # Check whether lead SNP is present in the extracted SNP list
    if ! grep -Fxq "${lead_SNP}" "${prefix}.snplist"; then
        echo "Warning: lead SNP ${lead_SNP} not found in ${prefix}.snplist, skipping."
        continue
    fi

    # Run PLINK LD
    "${PLINK_BIN}" \
        --bfile "${BFILE}" \
        --extract "${prefix}.snplist" \
        --ld-snp "${lead_SNP}" \
        --r2 \
        --ld-window 999999 \
        --ld-window-kb 1000000 \
        --ld-window-r2 0 \
        --out "${prefix}"

    echo "Finished: ${prefix}"
done

echo "End time: $(date)"
echo "All done."