#!/bin/bash
#SBATCH --job-name=plink_check_V2
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/plink-check-V2-%j.log

set -euo pipefail

module load plink-uoneasy/1.9b_6.21-x86_64


########################################
# 1. Input and output paths
########################################

PLINK_DIR=/gpfs01/home/shxyt1/data/sticklebacks/PLINK_sex_corrected

PLINK_MASTER=${PLINK_DIR}/stickleback_allchr_sex_corrected

GSEX_FILE=/gpfs01/home/shxyt1/data/sticklebacks/VCF/ploidy/Gsex_bcftools.tsv

X_BIM=${PLINK_DIR}/stickleback_NC_053230.1_sex_corrected.bim

Y_BIM=${PLINK_DIR}/stickleback_NC_053233.1_sex_corrected.bim

OUTDIR=/gpfs01/home/shxyt1/data/sticklebacks/plink_check_sex_corrected

mkdir -p "${OUTDIR}"


########################################
# 2. Check required files
########################################

for ext in bed bim fam; do
    if [[ ! -s "${PLINK_MASTER}.${ext}" ]]; then
        echo "ERROR: missing input file:"
        echo "${PLINK_MASTER}.${ext}"
        exit 1
    fi
done

for file in "${GSEX_FILE}" "${X_BIM}" "${Y_BIM}"; do
    if [[ ! -s "${file}" ]]; then
        echo "ERROR: missing input file:"
        echo "${file}"
        exit 1
    fi
done


########################################
# 3. Convert sex file for PLINK
#
# GSEX_FILE:
# QTL22002    F
# QTL22005    M
#
# PLINK update-sex:
# FID IID SEX
# SEX: 1=male, 2=female
########################################

PLINK_SEX_FILE=${OUTDIR}/sample_sex_for_plink.txt

awk '
BEGIN {
    OFS="\t"
}

$2=="M" {
    print $1, $1, 1
}

$2=="F" {
    print $1, $1, 2
}
' "${GSEX_FILE}" > "${PLINK_SEX_FILE}"

SEX_TOTAL=$(wc -l < "${PLINK_SEX_FILE}")
SEX_MALES=$(awk '$3==1 {n++} END{print n+0}' "${PLINK_SEX_FILE}")
SEX_FEMALES=$(awk '$3==2 {n++} END{print n+0}' "${PLINK_SEX_FILE}")

echo "Sex file samples: ${SEX_TOTAL}"
echo "Male samples: ${SEX_MALES}"
echo "Female samples: ${SEX_FEMALES}"

if [[ "${SEX_TOTAL}" -ne 352 ]]; then
    echo "ERROR: expected 352 samples in sex file"
    exit 1
fi

if [[ "${SEX_MALES}" -ne 103 ]]; then
    echo "ERROR: expected 103 males"
    exit 1
fi

if [[ "${SEX_FEMALES}" -ne 249 ]]; then
    echo "ERROR: expected 249 females"
    exit 1
fi


########################################
# 4. Add sex information to merged PLINK
########################################

STEP0_FULL=${OUTDIR}/step0_allchr_sex_updated

echo "[1] Updating sample sex in FAM file..."

plink \
    --bfile "${PLINK_MASTER}" \
    --update-sex "${PLINK_SEX_FILE}" \
    --allow-extra-chr \
    --keep-allele-order \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --memory 60000 \
    --make-bed \
    --out "${STEP0_FULL}"


########################################
# 5. Confirm sex was updated correctly
########################################

FAM_TOTAL=$(wc -l < "${STEP0_FULL}.fam")
FAM_MALES=$(awk '$5==1 {n++} END{print n+0}' "${STEP0_FULL}.fam")
FAM_FEMALES=$(awk '$5==2 {n++} END{print n+0}' "${STEP0_FULL}.fam")
FAM_UNKNOWN=$(awk '$5==0 {n++} END{print n+0}' "${STEP0_FULL}.fam")

echo "Samples in updated FAM: ${FAM_TOTAL}"
echo "Males in updated FAM: ${FAM_MALES}"
echo "Females in updated FAM: ${FAM_FEMALES}"
echo "Unknown sex in updated FAM: ${FAM_UNKNOWN}"

if [[ "${FAM_TOTAL}" -ne 352 ||
      "${FAM_MALES}" -ne 103 ||
      "${FAM_FEMALES}" -ne 249 ||
      "${FAM_UNKNOWN}" -ne 0 ]]; then

    echo "ERROR: sex information was not updated correctly"
    exit 1
fi


########################################
# 6. Create list of sex-chromosome SNPs
#
# NC_053230.1 = X/chr19
# NC_053233.1 = Y
########################################

SEX_CHR_SNPS=${OUTDIR}/sex_chromosome_variant_ids.txt

awk '{print $2}' \
    "${X_BIM}" \
    "${Y_BIM}" \
| sort -u \
> "${SEX_CHR_SNPS}"

SEX_CHR_COUNT=$(wc -l < "${SEX_CHR_SNPS}")

echo "Sex-chromosome variants to exclude from autosomal QC:"
echo "${SEX_CHR_COUNT}"

if [[ "${SEX_CHR_COUNT}" -eq 0 ]]; then
    echo "ERROR: sex chromosome SNP list is empty"
    exit 1
fi


########################################
# 7. Create autosomal dataset
########################################

STEP0_AUTOSOMES=${OUTDIR}/step0_autosomes

echo "[2] Creating autosomal dataset..."

plink \
    --bfile "${STEP0_FULL}" \
    --exclude "${SEX_CHR_SNPS}" \
    --allow-extra-chr \
    --keep-allele-order \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --memory 60000 \
    --make-bed \
    --out "${STEP0_AUTOSOMES}"


########################################
# 8. Basic QC on autosomes
#
# SNP missingness >5% removed
# Sample missingness >10% removed
########################################

STEP1_AUTOSOMAL_QC=${OUTDIR}/step1_autosomal_basicQC

echo "[3] Running autosomal missingness QC..."

plink \
    --bfile "${STEP0_AUTOSOMES}" \
    --geno 0.05 \
    --mind 0.10 \
    --allow-extra-chr \
    --keep-allele-order \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --memory 60000 \
    --make-bed \
    --out "${STEP1_AUTOSOMAL_QC}"


########################################
# 9. Create retained-sample list
########################################

KEEP_SAMPLES=${OUTDIR}/samples_passing_autosomal_QC.txt

awk 'BEGIN{OFS="\t"} {print $1, $2}' \
    "${STEP1_AUTOSOMAL_QC}.fam" \
    > "${KEEP_SAMPLES}"

RETAINED_SAMPLES=$(wc -l < "${KEEP_SAMPLES}")

echo "Samples retained after autosomal QC: ${RETAINED_SAMPLES}"


########################################
# 10. Missingness reports after QC
########################################

echo "[4] Generating missingness reports..."

plink \
    --bfile "${STEP1_AUTOSOMAL_QC}" \
    --missing \
    --allow-extra-chr \
    --out "${OUTDIR}/step2_missingness"


########################################
# 11. LD pruning for PCA and relatedness
#
# Common variants only for structure analyses
########################################

echo "[5] Running autosomal LD pruning..."

plink \
    --bfile "${STEP1_AUTOSOMAL_QC}" \
    --maf 0.05 \
    --indep-pairwise 50 5 0.2 \
    --allow-extra-chr \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --memory 60000 \
    --out "${OUTDIR}/step3_prunedSNPs"


########################################
# 12. PCA using pruned autosomal SNPs
########################################

echo "[6] Performing PCA..."

plink \
    --bfile "${STEP1_AUTOSOMAL_QC}" \
    --extract "${OUTDIR}/step3_prunedSNPs.prune.in" \
    --pca 20 \
    --allow-extra-chr \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --memory 60000 \
    --out "${OUTDIR}/step4_PCA"


########################################
# 13. Relatedness using pruned SNPs
########################################

echo "[7] Calculating relatedness..."

plink \
    --bfile "${STEP1_AUTOSOMAL_QC}" \
    --extract "${OUTDIR}/step3_prunedSNPs.prune.in" \
    --genome \
    --allow-extra-chr \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --memory 60000 \
    --out "${OUTDIR}/step5_relatedness"


########################################
# 14. Final autosomal GWAS dataset
########################################

FINAL_AUTOSOMES=${OUTDIR}/final_GWAS_ready_autosomes

echo "[8] Creating autosomal GWAS-ready dataset..."

plink \
    --bfile "${STEP1_AUTOSOMAL_QC}" \
    --allow-extra-chr \
    --keep-allele-order \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --memory 60000 \
    --make-bed \
    --out "${FINAL_AUTOSOMES}"


########################################
# 15. Final all-chromosome dataset
#
# Only samples failing autosomal QC are removed.
# Do not apply --geno globally because female Y
# ploidy=0 creates expected structural missingness.
########################################

FINAL_ALLCHR=${OUTDIR}/final_GWAS_ready_allchr

echo "[9] Creating all-chromosome sample-cleaned dataset..."

plink \
    --bfile "${STEP0_FULL}" \
    --keep "${KEEP_SAMPLES}" \
    --allow-extra-chr \
    --keep-allele-order \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --memory 60000 \
    --make-bed \
    --out "${FINAL_ALLCHR}"


########################################
# 16. Optional PLINK sex check
#
# --check-sex only works correctly when the
# sex chromosome is coded as X in the BIM.
########################################

X_CODED_VARIANTS=$(awk '$1=="X" {n++} END{print n+0}' "${STEP0_FULL}.bim")

if [[ "${X_CODED_VARIANTS}" -gt 0 ]]; then

    echo "[10] X chromosome detected. Running PLINK sex check..."

    plink \
        --bfile "${FINAL_ALLCHR}" \
        --check-sex \
        --allow-extra-chr \
        --out "${OUTDIR}/step6_sexcheck"

else

    echo "[10] Sex check skipped."
    echo "The sex chromosome is not coded as X in the merged BIM."
    echo "This does not affect the sex-aware SNP calling already completed."

fi


########################################
# 17. Final summary
########################################

AUTOSOMAL_VARIANTS=$(wc -l < "${FINAL_AUTOSOMES}.bim")
ALLCHR_VARIANTS=$(wc -l < "${FINAL_ALLCHR}.bim")
FINAL_SAMPLES=$(wc -l < "${FINAL_ALLCHR}.fam")

echo "========================================"
echo "PLINK QC pipeline completed successfully"
echo "========================================"
echo "Final samples: ${FINAL_SAMPLES}"
echo "Autosomal variants after QC: ${AUTOSOMAL_VARIANTS}"
echo "All-chromosome variants: ${ALLCHR_VARIANTS}"
echo
echo "Autosomal GWAS dataset:"
echo "${FINAL_AUTOSOMES}.{bed,bim,fam}"
echo
echo "All-chromosome sample-cleaned dataset:"
echo "${FINAL_ALLCHR}.{bed,bim,fam}"
echo
echo "PCA:"
echo "${OUTDIR}/step4_PCA.eigenvec"
echo
echo "Relatedness:"
echo "${OUTDIR}/step5_relatedness.genome"
echo "========================================"