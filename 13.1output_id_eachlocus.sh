#!/bin/bash
#SBATCH --job-name=.idperlocus
#SBATCH --partition=defq
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/idperlocus_%j.log

bfile=/gpfs01/home/shxyt1/data/sticklebacks/PLINK/stickleback_all_clean
in=/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_nine_main/locus_tbl2_for_ld.tsv
outdir=/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_nine_main/LD

mkdir -p "$outdir"

# Skip header, then read each locus line
tail -n +2 "$in" | while IFS=$'\t' read -r Element CHR locus_id lead_SNP window_start window_end; do
  outprefix="${outdir}/${Element}_chr${CHR}_locus${locus_id}_lead${lead_SNP}"

  plink --bfile "$bfile" \
    --chr "$CHR" \
    --from-bp "$window_start" \
    --to-bp "$window_end" \
    --r2 \
    --ld-snp "$lead_SNP" \
    --out "$outprefix"
done
