#!/bin/bash
#SBATCH --job-name=stickleback_simpleM
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2-00:00:00
#SBATCH --output=/gpfs01/home/shxyt1/app/simpleM/simpleM_%j.out
#SBATCH --error=/gpfs01/home/shxyt1/app/simpleM/simpleM_%j.err

set -euo pipefail


# ==========================================================
# Software
# ==========================================================

module load plink-uoneasy/1.9b_6.21-x86_64


# ==========================================================
# Input
# ==========================================================

BFILE="/gpfs01/home/shxyt1/data/sticklebacks/PLINK/stickleback_all_clean"


# ==========================================================
# Output
# ==========================================================

OUTDIR="/gpfs01/home/shxyt1/app/simpleM"

mkdir -p "${OUTDIR}"

TMPDIR="${OUTDIR}/tmp_${SLURM_JOB_ID}"

mkdir -p "${TMPDIR}"


# ==========================================================
# SimpleM parameters
# ==========================================================

BLOCK_SIZE=133

PCA_CUTOFF=0.995


# ==========================================================
# Chromosomes
# ==========================================================

CHROMS="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22"


# ==========================================================
# Output files
# ==========================================================

CHR_SUMMARY="${OUTDIR}/simpleM_chr_summary.txt"

BLOCK_SUMMARY="${OUTDIR}/simpleM_block_details.txt"

FINAL="${OUTDIR}/simpleM_final.txt"


echo -e \
"CHR\tN_SNP\tN_USED\tN_BLOCKS\tMeff\tMeff_ratio" \
> "${CHR_SUMMARY}"


echo -e \
"CHR\tBLOCK\tSTART_SNP\tEND_SNP\tSTART_BP\tEND_BP\tN_SNP\tN_USED\tMeff" \
> "${BLOCK_SUMMARY}"


# ==========================================================
# Check Python + numpy
# ==========================================================

if ! command -v python >/dev/null 2>&1
then

    echo "ERROR: python not found"

    exit 1

fi


python - <<'PY'
import numpy
print("Python/numpy check passed")
print("numpy version:", numpy.__version__)
PY


# ==========================================================
# Create temporary Python SimpleM calculator
# ==========================================================

PYTHON_SCRIPT="${TMPDIR}/simpleM_calc.py"


cat > "${PYTHON_SCRIPT}" <<'PYTHON'
#!/usr/bin/env python

import sys
import math
import numpy as np


# ==========================================================
# Arguments
# ==========================================================

traw_file = sys.argv[1]

chromosome = sys.argv[2]

block_size = int(sys.argv[3])

pca_cutoff = float(sys.argv[4])

block_output = sys.argv[5]

chr_output = sys.argv[6]


# ==========================================================
# SimpleM calculation for one SNP block
# ==========================================================

def calculate_meff(block):

    """
    block:
        list of tuples:
        (SNP_ID, BP, genotype_vector)

    Genotype vector is additive dosage 0/1/2.
    """

    n_original = len(block)

    if n_original == 0:
        return None


    # SNP x individual
    geno = np.asarray(
        [x[2] for x in block],
        dtype=float
    )


    # transpose:
    # individual x SNP
    geno = geno.T


    # ======================================================
    # Mean-impute missing genotype calls SNP-by-SNP
    # ======================================================

    for j in range(geno.shape[1]):

        col = geno[:, j]

        missing = np.isnan(col)

        if np.all(missing):

            continue

        if np.any(missing):

            mean_value = np.nanmean(col)

            col[missing] = mean_value

            geno[:, j] = col


    # ======================================================
    # Remove unusable / monomorphic SNPs
    # ======================================================

    keep = []

    for j in range(geno.shape[1]):

        col = geno[:, j]

        if np.any(np.isnan(col)):
            continue

        if np.std(col, ddof=1) <= 0:
            continue

        keep.append(j)


    if len(keep) == 0:

        return {
            "n_original": n_original,
            "n_used": 0,
            "meff": 0
        }


    geno = geno[:, keep]

    n_used = geno.shape[1]


    # ======================================================
    # One SNP
    # ======================================================

    if n_used == 1:

        return {
            "n_original": n_original,
            "n_used": 1,
            "meff": 1
        }


    # ======================================================
    # SNP-SNP Pearson correlation matrix
    # ======================================================

    corr = np.corrcoef(geno, rowvar=False)


    # Numerical clean-up
    corr = np.nan_to_num(
        corr,
        nan=0.0,
        posinf=0.0,
        neginf=0.0
    )

    np.fill_diagonal(corr, 1.0)


    # ======================================================
    # Eigenvalues
    #
    # eigvalsh is used because the correlation matrix
    # is symmetric.
    # ======================================================

    eigenvalues = np.linalg.eigvalsh(corr)


    # Match the SimpleM implementation:
    # absolute eigenvalues and descending order
    eigenvalues = np.sort(
        np.abs(eigenvalues)
    )[::-1]


    total = eigenvalues.sum()

    target = pca_cutoff * total

    cumulative = np.cumsum(eigenvalues)


    # smallest number of PCs explaining >= 99.5%
    meff = int(
        np.searchsorted(
            cumulative,
            target,
            side="left"
        ) + 1
    )


    return {
        "n_original": n_original,
        "n_used": n_used,
        "meff": meff
    }


# ==========================================================
# Read PLINK .traw file block-by-block
# ==========================================================

total_snps = 0

total_used = 0

total_meff = 0

n_blocks = 0


with open(traw_file, "r") as f:

    header = f.readline().strip().split()

    n_samples = len(header) - 6

    print(
        f"chr{chromosome}: {n_samples} individuals",
        flush=True
    )


    block = []


    for line in f:

        fields = line.strip().split()

        if len(fields) < 7:
            continue


        # PLINK .traw columns:
        #
        # CHR SNP CM POS COUNTED ALT sample1 sample2 ...
        #
        snp_id = fields[1]

        bp = int(fields[3])

        genotype_strings = fields[6:]


        genotype = []

        for value in genotype_strings:

            if value in ("NA", "NaN", "nan", "."):

                genotype.append(np.nan)

            else:

                genotype.append(float(value))


        block.append(
            (
                snp_id,
                bp,
                genotype
            )
        )


        # ==================================================
        # Process full block
        # ==================================================

        if len(block) == block_size:

            n_blocks += 1

            result = calculate_meff(block)


            start_snp = block[0][0]
            end_snp = block[-1][0]

            start_bp = block[0][1]
            end_bp = block[-1][1]


            with open(block_output, "a") as out:

                out.write(
                    f"{chromosome}\t"
                    f"{n_blocks}\t"
                    f"{start_snp}\t"
                    f"{end_snp}\t"
                    f"{start_bp}\t"
                    f"{end_bp}\t"
                    f"{result['n_original']}\t"
                    f"{result['n_used']}\t"
                    f"{result['meff']}\n"
                )


            total_snps += result["n_original"]

            total_used += result["n_used"]

            total_meff += result["meff"]


            if n_blocks % 500 == 0:

                print(
                    f"chr{chromosome}: "
                    f"{n_blocks} blocks processed; "
                    f"{total_snps} SNPs",
                    flush=True
                )


            block = []


    # ======================================================
    # Last incomplete block
    # ======================================================

    if len(block) > 0:

        n_blocks += 1

        result = calculate_meff(block)


        start_snp = block[0][0]
        end_snp = block[-1][0]

        start_bp = block[0][1]
        end_bp = block[-1][1]


        with open(block_output, "a") as out:

            out.write(
                f"{chromosome}\t"
                f"{n_blocks}\t"
                f"{start_snp}\t"
                f"{end_snp}\t"
                f"{start_bp}\t"
                f"{end_bp}\t"
                f"{result['n_original']}\t"
                f"{result['n_used']}\t"
                f"{result['meff']}\n"
            )


        total_snps += result["n_original"]

        total_used += result["n_used"]

        total_meff += result["meff"]


# ==========================================================
# Chromosome summary
# ==========================================================

ratio = (
    total_meff / total_used
    if total_used > 0
    else float("nan")
)


with open(chr_output, "a") as out:

    out.write(
        f"{chromosome}\t"
        f"{total_snps}\t"
        f"{total_used}\t"
        f"{n_blocks}\t"
        f"{total_meff}\t"
        f"{ratio:.8f}\n"
    )


print(
    f"chr{chromosome} finished: "
    f"N_SNP={total_snps}, "
    f"N_USED={total_used}, "
    f"blocks={n_blocks}, "
    f"Meff={total_meff}, "
    f"ratio={ratio:.4f}",
    flush=True
)

PYTHON


# ==========================================================
# Main chromosome loop
# ==========================================================

echo "========================================"
echo "Stickleback SimpleM"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Start: $(date)"
echo "BLOCK_SIZE: ${BLOCK_SIZE}"
echo "PCA_CUTOFF: ${PCA_CUTOFF}"
echo "========================================"


for CHR in ${CHROMS}
do

    echo
    echo "########################################"
    echo "Chromosome ${CHR}"
    echo "Start: $(date)"
    echo "########################################"


    PREFIX="${TMPDIR}/chr${CHR}"


    # ======================================================
    # Export additive genotype matrix
    #
    # .traw:
    # one SNP per row
    # one individual per genotype column
    # ======================================================

    plink \
        --bfile "${BFILE}" \
        --chr "${CHR}" \
        --recode A-transpose \
        --out "${PREFIX}"


    TRAW="${PREFIX}.traw"


    if [[ ! -s "${TRAW}" ]]
    then

        echo "ERROR: ${TRAW} not produced"

        exit 1

    fi


    # ======================================================
    # SimpleM
    # ======================================================

    python "${PYTHON_SCRIPT}" \
        "${TRAW}" \
        "${CHR}" \
        "${BLOCK_SIZE}" \
        "${PCA_CUTOFF}" \
        "${BLOCK_SUMMARY}" \
        "${CHR_SUMMARY}"


    # ======================================================
    # Delete temporary genotype export
    # ======================================================

    rm -f \
        "${PREFIX}.traw" \
        "${PREFIX}.log" \
        "${PREFIX}.nosex"


    echo "End chr${CHR}: $(date)"

done


# ==========================================================
# Genome-wide summary
# ==========================================================

TOTAL_SNP=$(awk '
    NR>1 {sum+=$2}
    END {print sum}
' "${CHR_SUMMARY}")


TOTAL_USED=$(awk '
    NR>1 {sum+=$3}
    END {print sum}
' "${CHR_SUMMARY}")


TOTAL_BLOCKS=$(awk '
    NR>1 {sum+=$4}
    END {print sum}
' "${CHR_SUMMARY}")


TOTAL_MEFF=$(awk '
    NR>1 {sum+=$5}
    END {print sum}
' "${CHR_SUMMARY}")


MEFF_RATIO=$(awk \
    -v m="${TOTAL_MEFF}" \
    -v n="${TOTAL_USED}" \
    'BEGIN {
        printf "%.8f", m/n
    }')


SIGNIFICANT=$(awk \
    -v m="${TOTAL_MEFF}" \
    'BEGIN {
        printf "%.12g", 0.05/m
    }')


SUGGESTIVE=$(awk \
    -v m="${TOTAL_MEFF}" \
    'BEGIN {
        printf "%.12g", 1/m
    }')


# ==========================================================
# Final output
# ==========================================================

{
    echo "========================================"
    echo "FINAL SimpleM SUMMARY"
    echo "========================================"
    echo
    echo "Block size            : ${BLOCK_SIZE} SNPs"
    echo "PCA cutoff            : ${PCA_CUTOFF}"
    echo
    echo "Total SNPs            : ${TOTAL_SNP}"
    echo "Expected SNPs         : 4549783"
    echo "SNPs used             : ${TOTAL_USED}"
    echo "Total blocks          : ${TOTAL_BLOCKS}"
    echo
    echo "Total Meff            : ${TOTAL_MEFF}"
    echo "Meff / SNP            : ${MEFF_RATIO}"
    echo
    echo "Suggestive 1/Meff     : ${SUGGESTIVE}"
    echo "Significant 0.05/Meff : ${SIGNIFICANT}"
    echo

    if [[ "${TOTAL_SNP}" -eq 4549783 ]]
    then
        echo "SNP CHECK: PASS"
    else
        echo "SNP CHECK: FAILED"
    fi

    echo
    echo "Finished: $(date)"
    echo "========================================"

} | tee "${FINAL}"


# ==========================================================
# Remove temporary directory
# ==========================================================

rm -rf "${TMPDIR}"