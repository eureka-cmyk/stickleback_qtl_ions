#!/bin/bash
#SBATCH --job-name=run_LD_P_V
#SBATCH --partition=defq
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/run_LD_P_V_%j.log
set -euo pipefail

module load bcftools-uoneasy/1.19-GCC-13.2.0
module load plink-uoneasy/1.9b_6.21-x86_64

PLINK="plink"

BFILE="/gpfs01/home/shxyt1/data/sticklebacks/PLINK/stickleback_all_clean"

LOCUS_FILE="/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_nine_main/six_main/regional plot/locus_tbl2_for_ld.tsv"

OUT_DIR="/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_nine_main/six_main"

mkdir -p "${OUT_DIR}"

if [ ! -f "${LOCUS_FILE}" ]; then
  echo "ERROR: locus file not found:"
  echo "${LOCUS_FILE}"
  exit 1
fi

if [ ! -f "${BFILE}.bed" ] || [ ! -f "${BFILE}.bim" ] || [ ! -f "${BFILE}.fam" ]; then
  echo "ERROR: PLINK bfile not found:"
  echo "${BFILE}"
  exit 1
fi

echo "Checking locus window sizes:"
awk -F '\t' '
NR==1 {
  for (i=1; i<=NF; i++) col[$i]=i
  print "Element\tCHR\tlocus_id\tlead_SNP\twindow_start\twindow_end\twindow_size"
  next
}
{
  size=$col["window_end"]-$col["window_start"]
  print $col["Element"]"\t"$col["CHR"]"\t"$col["locus_id"]"\t"$col["lead_SNP"]"\t"$col["window_start"]"\t"$col["window_end"]"\t"size
}
' "${LOCUS_FILE}"

echo "============================================"
echo "Start running LD for all loci from 150kb locus table"
echo "============================================"

awk -F '\t' '
NR==1 {
  for (i=1; i<=NF; i++) col[$i]=i
  next
}
{
  print $col["Element"], $col["CHR"], $col["locus_id"], $col["lead_SNP"], $col["window_start"], $col["window_end"]
}
' "${LOCUS_FILE}" | while read -r ELEMENT CHR LOCUS_ID LEAD_SNP START_BP END_BP
do
  echo "============================================"
  echo "Running LD: ${ELEMENT} chr${CHR} locus${LOCUS_ID} lead ${LEAD_SNP}"
  echo "Region: ${START_BP}-${END_BP}"

  OUT_PREFIX="${OUT_DIR}/${ELEMENT}_chr${CHR}_locus${LOCUS_ID}_lead${LEAD_SNP}"

  ${PLINK} \
    --bfile "${BFILE}" \
    --chr "${CHR}" \
    --from-bp "${START_BP}" \
    --to-bp "${END_BP}" \
    --ld-snp "${LEAD_SNP}" \
    --r2 \
    --ld-window 999999 \
    --ld-window-kb 1000000 \
    --ld-window-r2 0 \
    --out "${OUT_PREFIX}"

  echo "Saved: ${OUT_PREFIX}.ld"

done

echo "============================================"
echo "Done. LD files saved in:"
echo "${OUT_DIR}"