#!/bin/bash

set -uo pipefail

WORK="/gpfs01/home/shxyt1/app/GEC/stickleback_Me_run"

BIM="/gpfs01/home/shxyt1/data/sticklebacks/PLINK/stickleback_all_clean.bim"

WHOLE_DIR="${WORK}/by_chromosome"
REGION_DIR="${WORK}/by_region"
STATUS="${REGION_DIR}/region_status.txt"

OUT="${WORK}/final_GEC_summary.txt"

# ==========================================================
# Whole chromosomes that successfully ran without splitting
# ==========================================================

WHOLE_CHRS="1 3 5 7 8 9 11 12 13 14 15 16 18 19 20 22"


# ==========================================================
# Extract Me from a KGGSEE log
# ==========================================================

get_me () {

    LOG="$1"

    awk '
        /The effective number of independent tests for all the/ {
            for(i=1;i<=NF;i++){
                if($i=="is"){
                    print $(i+1)
                    exit
                }
            }
        }
    ' "${LOG}"
}


# ==========================================================
# Output header
# ==========================================================

echo -e "TYPE\tCHR\tSTART\tEND\tN_SNP\tMe" > "${OUT}"

TOTAL_ME=0
TOTAL_SNP=0
N_WHOLE=0
N_REGION=0


# ==========================================================
# 1. Whole chromosomes
# ==========================================================

for CHR in ${WHOLE_CHRS}
do

    LOG="${WHOLE_DIR}/chr${CHR}.log"

    if [[ ! -s "${LOG}" ]]; then
        echo "ERROR: Missing whole chromosome log:"
        echo "${LOG}"
        exit 1
    fi

    ME=$(get_me "${LOG}")

    if [[ -z "${ME}" ]]; then
        echo "ERROR: Cannot extract Me from:"
        echo "${LOG}"
        exit 1
    fi

    N_SNP=$(awk -v c="${CHR}" '$1==c {n++} END{print n+0}' "${BIM}")

    echo -e "WHOLE\t${CHR}\t1\tNA\t${N_SNP}\t${ME}" >> "${OUT}"

    TOTAL_ME=$(awk -v a="${TOTAL_ME}" -v b="${ME}" 'BEGIN{printf "%.10f",a+b}')
    TOTAL_SNP=$((TOTAL_SNP + N_SNP))
    N_WHOLE=$((N_WHOLE + 1))

done


# ==========================================================
# Function: identify the four obsolete regions
# ==========================================================

is_obsolete () {

    CHR="$1"
    START="$2"
    END="$3"

    [[ "${CHR}" == "4"  && "${START}" == "18000001" && "${END}" == "20000000" ]] && return 0
    [[ "${CHR}" == "4"  && "${START}" == "20000001" && "${END}" == "22000000" ]] && return 0

    [[ "${CHR}" == "21" && "${START}" == "8000001"  && "${END}" == "10000000" ]] && return 0
    [[ "${CHR}" == "21" && "${START}" == "10000001" && "${END}" == "12000000" ]] && return 0

    return 1
}


# ==========================================================
# 2. Original regional results
#    Skip the four regions replaced by optimized boundaries
# ==========================================================

while read CHR START END N_SNP STATE
do

    [[ "${CHR}" == "CHR" ]] && continue
    [[ "${STATE}" != "SUCCESS" ]] && continue

    if is_obsolete "${CHR}" "${START}" "${END}"; then

        echo "Skipping obsolete region: chr${CHR}:${START}-${END}"
        continue

    fi

    NAME="chr${CHR}_${START}_${END}"

    LOG="${REGION_DIR}/${NAME}/${NAME}.log"

    if [[ ! -s "${LOG}" ]]; then
        echo "ERROR: Missing regional log:"
        echo "${LOG}"
        exit 1
    fi

    ME=$(get_me "${LOG}")

    if [[ -z "${ME}" ]]; then
        echo "ERROR: Cannot extract Me from:"
        echo "${LOG}"
        exit 1
    fi

    echo -e "REGION\t${CHR}\t${START}\t${END}\t${N_SNP}\t${ME}" >> "${OUT}"

    TOTAL_ME=$(awk -v a="${TOTAL_ME}" -v b="${ME}" 'BEGIN{printf "%.10f",a+b}')
    TOTAL_SNP=$((TOTAL_SNP + N_SNP))
    N_REGION=$((N_REGION + 1))

done < "${STATUS}"


# ==========================================================
# 3. Add the four optimized regions
# ==========================================================

OPTIMIZED=(
"4 18000001 20500000"
"4 20500001 22000000"
"21 8000001 9900000"
"21 9900001 12000000"
)


for R in "${OPTIMIZED[@]}"
do

    read CHR START END <<< "${R}"

    NAME="chr${CHR}_${START}_${END}"

    LOG="${REGION_DIR}/${NAME}/${NAME}.log"

    if [[ ! -s "${LOG}" ]]; then
        echo "ERROR: Missing optimized log:"
        echo "${LOG}"
        exit 1
    fi

    ME=$(get_me "${LOG}")

    if [[ -z "${ME}" ]]; then
        echo "ERROR: Cannot extract Me from:"
        echo "${LOG}"
        exit 1
    fi

    N_SNP=$(awk \
        -v c="${CHR}" \
        -v s="${START}" \
        -v e="${END}" \
        '$1==c && $4>=s && $4<=e {n++}
         END{print n+0}' \
        "${BIM}")

    echo -e "OPTIMIZED\t${CHR}\t${START}\t${END}\t${N_SNP}\t${ME}" >> "${OUT}"

    TOTAL_ME=$(awk -v a="${TOTAL_ME}" -v b="${ME}" 'BEGIN{printf "%.10f",a+b}')
    TOTAL_SNP=$((TOTAL_SNP + N_SNP))
    N_REGION=$((N_REGION + 1))

done


# ==========================================================
# Final thresholds
# ==========================================================

SIGNIFICANT=$(awk -v me="${TOTAL_ME}" 'BEGIN{printf "%.12g",0.05/me}')
SUGGESTIVE=$(awk -v me="${TOTAL_ME}" 'BEGIN{printf "%.12g",1/me}')


echo
echo "========================================"
echo "FINAL GEC SUMMARY"
echo "========================================"

echo "Whole chromosomes : ${N_WHOLE}"
echo "Regional blocks   : ${N_REGION}"
echo
echo "Total SNPs        : ${TOTAL_SNP}"
echo "Expected SNPs     : 4549783"
echo
echo "Total Me          : ${TOTAL_ME}"
echo
echo "Suggestive 1/Me   : ${SUGGESTIVE}"
echo "Significant .05/Me: ${SIGNIFICANT}"
echo


if [[ "${TOTAL_SNP}" -eq 4549783 ]]; then

    echo "SNP CHECK: PASS"

else

    echo "SNP CHECK: FAILED"

fi


echo
echo "Detailed results:"
echo "${OUT}"
echo "========================================"