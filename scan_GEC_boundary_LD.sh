#!/bin/bash
#SBATCH --job-name=scan_GEC_boundary
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=2-00:00:00
#SBATCH --output=/gpfs01/home/shxyt1/app/GEC/LD_boundary_scan/scan_boundary_%j.out
#SBATCH --error=/gpfs01/home/shxyt1/app/GEC/LD_boundary_scan/scan_boundary_%j.err

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

OUTDIR="/gpfs01/home/shxyt1/app/GEC/LD_boundary_scan"

mkdir -p "${OUTDIR}"

ALL_RESULT="${OUTDIR}/all_candidate_boundaries.txt"
BEST_RESULT="${OUTDIR}/best_candidate_boundaries.txt"


# ==========================================================
# Settings
# ==========================================================

CHROMS="2 4 6 10 17 21"

# Original GEC region size
BLOCK_SIZE=2000000

# Search +/- 500 kb around original boundary
SCAN_RANGE=500000

# Candidate boundary every 50 kb
STEP=50000

# For each candidate, examine SNPs within +/- 500 kb
FLANK=500000

# Closest SNPs from each side
N_SIDE=500


# ==========================================================
# Output headers
# ==========================================================

echo -e \
"CHR\tOriginal_Boundary_BP\tCandidate_BP\tShift_BP\tLeft_N\tRight_N\tN_cross_pairs\tMean_R2\tMax_R2\tProp_R2_gt_0.1\tProp_R2_gt_0.2" \
> "${ALL_RESULT}"

echo -e \
"CHR\tOriginal_Boundary_BP\tBest_Boundary_BP\tShift_BP\tLeft_N\tRight_N\tN_cross_pairs\tMean_R2\tMax_R2\tProp_R2_gt_0.1\tProp_R2_gt_0.2" \
> "${BEST_RESULT}"


echo "========================================"
echo "Scanning GEC artificial boundaries"
echo "Start: $(date)"
echo "Job ID: ${SLURM_JOB_ID:-NA}"
echo "Search range: +/- ${SCAN_RANGE} bp"
echo "Step: ${STEP} bp"
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

    CHR_SNPS="${CHRDIR}/chr${CHR}.position_snps.txt"


    # ------------------------------------------------------
    # Position + SNP ID
    # ------------------------------------------------------

    awk -v c="${CHR}" '
        $1==c {
            print $4 "\t" $2
        }
    ' "${BIM}" | sort -n -k1,1 > "${CHR_SNPS}"


    MAX_BP=$(tail -n 1 "${CHR_SNPS}" | awk '{print $1}')

    echo "Maximum BP: ${MAX_BP}"


    # ======================================================
    # Original artificial boundaries
    # ======================================================

    ORIGINAL=${BLOCK_SIZE}

    while [[ ${ORIGINAL} -lt ${MAX_BP} ]]
    do

        echo
        echo "========================================"
        echo "Original boundary: chr${CHR}:${ORIGINAL}"
        echo "========================================"


        ORIGINAL_DIR="${CHRDIR}/boundary_${ORIGINAL}"

        mkdir -p "${ORIGINAL_DIR}"


        CANDIDATE_RESULT="${ORIGINAL_DIR}/candidate_results.txt"

        echo -e \
"CHR\tOriginal_Boundary_BP\tCandidate_BP\tShift_BP\tLeft_N\tRight_N\tN_cross_pairs\tMean_R2\tMax_R2\tProp_R2_gt_0.1\tProp_R2_gt_0.2" \
        > "${CANDIDATE_RESULT}"


        # --------------------------------------------------
        # Candidate search range
        # --------------------------------------------------

        SCAN_START=$((ORIGINAL - SCAN_RANGE))
        SCAN_END=$((ORIGINAL + SCAN_RANGE))


        if [[ ${SCAN_START} -lt 1 ]]; then
            SCAN_START=1
        fi


        if [[ ${SCAN_END} -ge ${MAX_BP} ]]; then
            SCAN_END=$((MAX_BP - 1))
        fi


        CANDIDATE=${SCAN_START}


        # ==================================================
        # Candidate boundary loop
        # ==================================================

        while [[ ${CANDIDATE} -le ${SCAN_END} ]]
        do

            SHIFT=$((CANDIDATE - ORIGINAL))


            echo
            echo "Candidate: chr${CHR}:${CANDIDATE}   shift=${SHIFT}"


            PREFIX="${ORIGINAL_DIR}/candidate_${CANDIDATE}"

            LEFT_LIST="${PREFIX}.left.txt"
            RIGHT_LIST="${PREFIX}.right.txt"
            EXTRACT_LIST="${PREFIX}.extract.txt"


            LEFT_START=$((CANDIDATE - FLANK))
            RIGHT_END=$((CANDIDATE + FLANK))


            if [[ ${LEFT_START} -lt 1 ]]; then
                LEFT_START=1
            fi


            if [[ ${RIGHT_END} -gt ${MAX_BP} ]]; then
                RIGHT_END=${MAX_BP}
            fi


            # ==============================================
            # Closest N_SIDE SNPs on LEFT
            # ==============================================

            awk \
                -v lo="${LEFT_START}" \
                -v b="${CANDIDATE}" \
                '$1>=lo && $1<=b {print}' \
                "${CHR_SNPS}" \
            | tail -n "${N_SIDE}" \
            | awk '{print $2}' \
            > "${LEFT_LIST}"


            # ==============================================
            # Closest N_SIDE SNPs on RIGHT
            # ==============================================

            awk \
                -v b="${CANDIDATE}" \
                -v hi="${RIGHT_END}" \
                '$1>b && $1<=hi {print}' \
                "${CHR_SNPS}" \
            | head -n "${N_SIDE}" \
            | awk '{print $2}' \
            > "${RIGHT_LIST}"


            LEFT_N=$(wc -l < "${LEFT_LIST}")
            RIGHT_N=$(wc -l < "${RIGHT_LIST}")


            echo "Left N : ${LEFT_N}"
            echo "Right N: ${RIGHT_N}"


            # ----------------------------------------------
            # No SNPs on one side
            # ----------------------------------------------

            if [[ ${LEFT_N} -eq 0 || ${RIGHT_N} -eq 0 ]]; then

                LINE="${CHR}\t${ORIGINAL}\t${CANDIDATE}\t${SHIFT}\t${LEFT_N}\t${RIGHT_N}\t0\tNA\tNA\tNA\tNA"

                echo -e "${LINE}" >> "${CANDIDATE_RESULT}"
                echo -e "${LINE}" >> "${ALL_RESULT}"

                rm -f "${LEFT_LIST}" "${RIGHT_LIST}"

                CANDIDATE=$((CANDIDATE + STEP))

                continue

            fi


            # ==============================================
            # Extract SNPs
            # ==============================================

            cat "${LEFT_LIST}" "${RIGHT_LIST}" \
            | sort -u \
            > "${EXTRACT_LIST}"


            # ==============================================
            # PLINK LD
            # ==============================================

            plink \
                --bfile "${BFILE}" \
                --extract "${EXTRACT_LIST}" \
                --r2 gz \
                --ld-window 999999 \
                --ld-window-kb 1000 \
                --ld-window-r2 0 \
                --threads "${SLURM_CPUS_PER_TASK:-4}" \
                --out "${PREFIX}" \
                >/dev/null 2>&1


            LD_FILE="${PREFIX}.ld.gz"


            if [[ ! -s "${LD_FILE}" ]]; then

                echo "WARNING: no LD file"

                LINE="${CHR}\t${ORIGINAL}\t${CANDIDATE}\t${SHIFT}\t${LEFT_N}\t${RIGHT_N}\t0\tNA\tNA\tNA\tNA"

                echo -e "${LINE}" >> "${CANDIDATE_RESULT}"
                echo -e "${LINE}" >> "${ALL_RESULT}"

                rm -f "${LEFT_LIST}" \
                      "${RIGHT_LIST}" \
                      "${EXTRACT_LIST}"

                CANDIDATE=$((CANDIDATE + STEP))

                continue

            fi


            # ==============================================
            # Cross-boundary LD only
            # ==============================================

            STATS=$(zcat "${LD_FILE}" | \
            awk -v b="${CANDIDATE}" '

                NR>1 {

                    bp1=$2
                    bp2=$5
                    r2=$7

                    if ((bp1<=b && bp2>b) || (bp2<=b && bp1>b)) {

                        n++
                        sum+=r2

                        if (n==1 || r2>max) {
                            max=r2
                        }

                        if (r2>0.1) {
                            gt01++
                        }

                        if (r2>0.2) {
                            gt02++
                        }
                    }
                }

                END {

                    if (n>0) {

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


            LINE="${CHR}\t${ORIGINAL}\t${CANDIDATE}\t${SHIFT}\t${LEFT_N}\t${RIGHT_N}\t${STATS}"

            echo -e "${LINE}" >> "${CANDIDATE_RESULT}"
            echo -e "${LINE}" >> "${ALL_RESULT}"


            echo "Stats: ${STATS}"


            # ==============================================
            # Remove temporary files
            # ==============================================

            rm -f "${LD_FILE}"
            rm -f "${PREFIX}.log"
            rm -f "${PREFIX}.nosex"

            rm -f "${LEFT_LIST}"
            rm -f "${RIGHT_LIST}"
            rm -f "${EXTRACT_LIST}"


            CANDIDATE=$((CANDIDATE + STEP))

        done


        # ==================================================
        # Pick best candidate
        #
        # 1. Lowest mean R2
        # 2. If tied, lowest proportion R2 > 0.2
        # 3. If still tied, lowest proportion R2 > 0.1
        # ==================================================

        BEST=$(awk '

            NR==1 {
                next
            }

            $8!="NA" {

                if (!found) {

                    found=1

                    chr=$1
                    orig=$2
                    cand=$3
                    shift=$4
                    ln=$5
                    rn=$6
                    np=$7
                    mean=$8
                    max=$9
                    p01=$10
                    p02=$11

                } else if ($8 < mean || ($8 == mean && $11 < p02) || ($8 == mean && $11 == p02 && $10 < p01)) {

                    chr=$1
                    orig=$2
                    cand=$3
                    shift=$4
                    ln=$5
                    rn=$6
                    np=$7
                    mean=$8
                    max=$9
                    p01=$10
                    p02=$11
                }
            }

            END {

                if (found) {

                    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s",
                           chr,
                           orig,
                           cand,
                           shift,
                           ln,
                           rn,
                           np,
                           mean,
                           max,
                           p01,
                           p02
                }
            }

        ' "${CANDIDATE_RESULT}")


        if [[ -n "${BEST}" ]]; then

            echo -e "${BEST}" >> "${BEST_RESULT}"

            echo
            echo "BEST:"
            echo -e "${BEST}"

        else

            echo "WARNING: no valid candidate found"

        fi


        ORIGINAL=$((ORIGINAL + BLOCK_SIZE))

    done

done


echo
echo "========================================"
echo "Boundary scan finished"
echo "Finished: $(date)"
echo "========================================"

echo
echo "Best candidate boundaries:"
column -t "${BEST_RESULT}"