#!/bin/bash
#SBATCH --job-name=GEC_optimized
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs01/home/shxyt1/app/GEC/GEC_optimized_%j.out
#SBATCH --error=/gpfs01/home/shxyt1/app/GEC/GEC_optimized_%j.err

set -uo pipefail

echo "========================================"
echo "Optimized GEC regions"
echo "Job ID: ${SLURM_JOB_ID:-NA}"
echo "Start: $(date)"
echo "Host: $(hostname)"
echo "========================================"

module load java-uoneasy/17.0.6

GEC_DIR="/gpfs01/home/shxyt1/app/GEC"
JAR="${GEC_DIR}/kggsee.jar"

WORK="${GEC_DIR}/stickleback_Me_run"

BIM="/gpfs01/home/shxyt1/data/sticklebacks/PLINK/stickleback_all_clean.bim"

VCF_DIR="${WORK}/vcf"
PFILE_DIR="${WORK}/pfiles_by_chr"
OUTDIR="${WORK}/by_region"

mkdir -p "${OUTDIR}"


# ==========================================================
# Optimized regions
# ==========================================================

REGIONS=(
"4 18000001 20500000"
"4 20500001 22000000"
"21 8000001 9900000"
"21 9900001 12000000"
)


for REGION in "${REGIONS[@]}"
do

    read CHR START END <<< "${REGION}"

    REGION_NAME="chr${CHR}_${START}_${END}"
    REGION_DIR="${OUTDIR}/${REGION_NAME}"

    mkdir -p "${REGION_DIR}"

    VCF="${VCF_DIR}/stickleback_gwas.chr${CHR}.vcf.gz"
    PFILE="${PFILE_DIR}/chr${CHR}.pvalue.txt"

    OUT="${REGION_DIR}/${REGION_NAME}"

    N_SNP=$(awk \
        -v c="${CHR}" \
        -v s="${START}" \
        -v e="${END}" \
        '$1==c && $4>=s && $4<=e {n++}
         END {print n+0}' \
        "${BIM}")

    echo
    echo "========================================"
    echo "Region: chr${CHR}:${START}-${END}"
    echo "SNPs: ${N_SNP}"
    echo "Start: $(date)"
    echo "========================================"


    # ------------------------------------------------------
    # Skip if already successfully completed
    # ------------------------------------------------------

    if [[ -s "${OUT}.effective.size.txt.gz" ]] && \
       gzip -t "${OUT}.effective.size.txt.gz" 2>/dev/null
    then

        N_LINES=$(zcat "${OUT}.effective.size.txt.gz" | wc -l)

        if [[ "${N_LINES}" -gt 1 ]]; then
            echo "Already SUCCESS -> skipping"
            continue
        fi
    fi


    # ------------------------------------------------------
    # Clean incomplete output only
    # ------------------------------------------------------

    rm -f "${OUT}.effective.size.txt.gz"
    rm -f "${OUT}.log"
    rm -rf "${REGION_DIR}/VCFRefhg19"

    cd "${REGION_DIR}"


    # ------------------------------------------------------
    # GEC
    # ------------------------------------------------------

    java -Xmx48g -jar "${JAR}" \
        --var-gec \
        --nt "${SLURM_CPUS_PER_TASK:-8}" \
        --chrom 1-22 \
        --filter-maf-le 0 \
        --hwe-all 0 \
        --ld-block-r2 0.10 \
        --regions-in "chr${CHR}:${START}-${END}" \
        --p-value-cutoff 0.05 \
        --vcf-ref "${VCF}" \
        --pfile "${PFILE}" \
        --out "${OUT}"

    EXITCODE=$?


    # ------------------------------------------------------
    # Validate
    # ------------------------------------------------------

    if [[ ${EXITCODE} -eq 0 ]] && \
       [[ -s "${OUT}.effective.size.txt.gz" ]] && \
       gzip -t "${OUT}.effective.size.txt.gz" 2>/dev/null
    then

        N_LINES=$(zcat "${OUT}.effective.size.txt.gz" | wc -l)

        if [[ "${N_LINES}" -gt 1 ]]; then
            echo "Region SUCCESS"
        else
            echo "Region FAILED: empty result"
            rm -f "${OUT}.effective.size.txt.gz"
        fi

    else

        echo "Region FAILED"
        rm -f "${OUT}.effective.size.txt.gz"

    fi


    rm -rf "${REGION_DIR}/VCFRefhg19"

    echo "End: $(date)"

done


echo
echo "========================================"
echo "Optimized regions finished"
echo "Finished: $(date)"
echo "========================================"