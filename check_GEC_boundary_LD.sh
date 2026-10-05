#!/bin/bash
#SBATCH --job-name=GEC_boundary_LD
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs01/home/shxyt1/app/GEC/LD_boundary/boundary_LD_%j.out
#SBATCH --error=/gpfs01/home/shxyt1/app/GEC/LD_boundary/boundary_LD_%j.err

set -uo pipefail

# ==========================================================
# Software
# ==========================================================

module load plink-uoneasy/1.9b_6.21-x86_64


# ==========================================================
# Paths
# ==========================================================

BFILE="/gpfs01/home/shxyt1/data/sticklebacks/PLINK/stickleback_all_clean"
BIM="${BFILE}.bim"

OUTDIR="/gpfs01/home/shxyt1/app/GEC/LD_boundary"
mkdir -p "${OUTDIR}"

SUMMARY="${OUTDIR}/boundary_LD_summary.txt"


# ==========================================================
# Settings
# ==========================================================

# Chromosomes which were divided into 2-Mb GEC regions
CHROMS="2 4 6 10 17 21"

# Artificial GEC block size
BLOCK_SIZE=2000000

# Look up to 500 kb on either side of each boundary
FLANK=500000

# Use the closest 500 SNPs from each side
N_SIDE=500


# ==========================================================
# Summary header
# ==========================================================

echo -e \
"CHR\tBoundary_BP\tBoundary_Mb\tLeft_N\tRight_N\tN_cross_pairs\tMean_R2\tMax_R2\tProp_R2_gt_0.1\tProp_R2_gt_0.2" \
> "${SUMMARY}"


echo "========================================"
echo "GEC boundary LD analysis"
echo "Start: $(date)"
echo "Job ID: ${SLURM_JOB_ID:-NA}"
echo "========================================"


# ==========================================================
# Chromosome loop
# ==========================================================

for CHR in ${CHROMS}
do

    echo
    echo "########################################"
    echo "Chromosome ${CHR}"
    echo "########################################"

    CHRDIR="${OUTDIR}/chr${CHR}"
    mkdir -p "${CHRDIR}"

    # Make a simple position + SNP-ID file once for this chromosome
    CHR_SNPS="${CHRDIR}/chr${CHR}.position_snps.txt"

    awk -v c="${CHR}" '
        $1==c {
            print $4 "\t" $2
        }
    ' "${BIM}" | sort -n -k1,1 > "${CHR_SNPS}"


    MAX_BP=$(tail -n 1 "${CHR_SNPS}" | awk '{print $1}')

    echo "Maximum BP: ${MAX_BP}"

    # First artificial boundary = 2 Mb
    BOUNDARY=${BLOCK_SIZE}


    # ======================================================
    # Boundary loop
    # ======================================================

    while [[ ${BOUNDARY} -lt ${MAX_BP} ]]
    do

        BOUNDARY_MB=$((BOUNDARY / 1000000))

        echo
        echo "----------------------------------------"
        echo "chr${CHR} boundary: ${BOUNDARY_MB} Mb"
        echo "----------------------------------------"

        PREFIX="${CHRDIR}/chr${CHR}_boundary_${BOUNDARY}"

        LEFT_LIST="${PREFIX}.left.txt"
        RIGHT_LIST="${PREFIX}.right.txt"
        EXTRACT_LIST="${PREFIX}.extract.txt"

        LEFT_START=$((BOUNDARY - FLANK))

        if [[ ${LEFT_START} -lt 1 ]]; then
            LEFT_START=1
        fi

        RIGHT_END=$((BOUNDARY + FLANK))

        if [[ ${RIGHT_END} -gt ${MAX_BP} ]]; then
            RIGHT_END=${MAX_BP}
        fi


        # ==================================================
        # Closest N_SIDE SNPs immediately LEFT of boundary
        # ==================================================

        awk \
            -v lo="${LEFT_START}" \
            -v b="${BOUNDARY}" \
            '$1>=lo && $1<=b {print}' \
            "${CHR_SNPS}" \
        | tail -n "${N_SIDE}" \
        | awk '{print $2}' \
        > "${LEFT_LIST}"


        # ==================================================
        # Closest N_SIDE SNPs immediately RIGHT of boundary
        # ==================================================

        awk \
            -v b="${BOUNDARY}" \
            -v hi="${RIGHT_END}" \
            '$1>b && $1<=hi {print}' \
            "${CHR_SNPS}" \
        | head -n "${N_SIDE}" \
        | awk '{print $2}' \
        > "${RIGHT_LIST}"


        LEFT_N=$(wc -l < "${LEFT_LIST}")
        RIGHT_N=$(wc -l < "${RIGHT_LIST}")


        echo "Left SNPs : ${LEFT_N}"
        echo "Right SNPs: ${RIGHT_N}"


        if [[ ${LEFT_N} -eq 0 || ${RIGHT_N} -eq 0 ]]; then

            echo "WARNING: no SNPs on one side of boundary"

            echo -e \
"${CHR}\t${BOUNDARY}\t${BOUNDARY_MB}\t${LEFT_N}\t${RIGHT_N}\t0\tNA\tNA\tNA\tNA" \
            >> "${SUMMARY}"

            BOUNDARY=$((BOUNDARY + BLOCK_SIZE))

            continue

        fi


        cat "${LEFT_LIST}" "${RIGHT_LIST}" \
        | sort -u \
        > "${EXTRACT_LIST}"


        # ==================================================
        # Calculate pairwise LD
        #
        # 500 kb each side -> max possible cross-boundary
        # distance is <= 1 Mb.
        # ==================================================

        plink \
            --bfile "${BFILE}" \
            --extract "${EXTRACT_LIST}" \
            --r2 gz \
            --ld-window 999999 \
            --ld-window-kb 1000 \
            --ld-window-r2 0 \
            --threads "${SLURM_CPUS_PER_TASK:-4}" \
            --out "${PREFIX}"


        LD_FILE="${PREFIX}.ld.gz"


        if [[ ! -s "${LD_FILE}" ]]; then

            echo "ERROR: LD file not produced"

            echo -e \
"${CHR}\t${BOUNDARY}\t${BOUNDARY_MB}\t${LEFT_N}\t${RIGHT_N}\t0\tNA\tNA\tNA\tNA" \
            >> "${SUMMARY}"

            BOUNDARY=$((BOUNDARY + BLOCK_SIZE))

            continue

        fi


        # ==================================================
        # Keep ONLY pairs crossing the artificial boundary
        #
        # left SNP:  BP <= boundary
        # right SNP: BP >  boundary
        # ==================================================

        STATS=$(zcat "${LD_FILE}" | \
        awk -v b="${BOUNDARY}" '

            NR>1 {

                bp1=$2
                bp2=$5
                r2=$7

                # Pair must lie on opposite sides
                if ((bp1<=b && bp2>b) || (bp2<=b && bp1>b)) {

                    n++
                    sum+=r2

                    if(n==1 || r2>max){
                        max=r2
                    }

                    if(r2>0.1){
                        gt01++
                    }

                    if(r2>0.2){
                        gt02++
                    }
                }
            }

            END {

                if(n>0){

                    printf "%d\t%.8f\t%.8f\t%.8f\t%.8f",
                           n,
                           sum/n,
                           max,
                           gt01/n,
                           gt02/n

                } else {

                    printf "0\tNA\tNA\tNA\tNA"

                }
            }
        ')


        echo -e \
"${CHR}\t${BOUNDARY}\t${BOUNDARY_MB}\t${LEFT_N}\t${RIGHT_N}\t${STATS}" \
        >> "${SUMMARY}"


        echo "Cross-boundary statistics:"
        echo "${STATS}"


        # ==================================================
        # Remove large temporary LD-pair output.
        #
        # Summary statistics are retained.
        # PLINK log is retained for reproducibility.
        # ==================================================

        rm -f "${LD_FILE}"

        rm -f "${PREFIX}.nosex"

        rm -f "${LEFT_LIST}"
        rm -f "${RIGHT_LIST}"
        rm -f "${EXTRACT_LIST}"


        BOUNDARY=$((BOUNDARY + BLOCK_SIZE))

    done

done


echo
echo "========================================"
echo "Boundary LD analysis finished"
echo "Finished: $(date)"
echo "========================================"

echo
echo "Results:"
cat "${SUMMARY}"