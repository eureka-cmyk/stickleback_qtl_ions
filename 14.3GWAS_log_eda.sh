#!/bin/bash
#SBATCH --job-name=GWAS_EDA
#SBATCH --partition=defq
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --array=1-11
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/GWAS_EDA_%j.log

# PCS=/gpfs01/home/shxyt1/data/sticklebacks/pheno/traits/pcs5.txt
# EDA=/gpfs01/home/shxyt1/data/sticklebacks/pheno/traits/EDA_kmeans_1D_all_windows.txt
# OUT=/gpfs01/home/shxyt1/data/sticklebacks/pheno/traits/pc5_eda.txt

# awk 'BEGIN{FS=OFS="\t"}
# NR==FNR{
#   if(FNR==1){
#     for(i=1;i<=NF;i++) if($i=="EDA25_class") c=i
#     next
#   }
#   if(c>0 && $c!="" && $c!="NA") a[$1]=$c   # $1 是 IID
#   next
# }
# {
#   iid=$2  # pcs5.txt: FID IID ...
#   if(iid in a) print $0, a[iid]
# }' "$EDA" "$PCS" > "$OUT"


# head -n 3 "$OUT"
# awk -F'\t' 'NR==1{print "Columns in OUT =",NF}' "$OUT"



echo "Start time: $(date)"
echo "Hostname: $(hostname)"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"

# ============================================================
# 1. 接受 sbatch 输入变量（用户可在提交时覆盖）
# ============================================================

GCTA_BIN=${GCTA_BIN:-"/gpfs01/home/shxyt1/app/gcta/gcta-1.95.0-linux-kernel-3-x86_64/squashfs-root/usr/bin/gcta64"}
BFILE=${BFILE:-"/gpfs01/home/shxyt1/data/sticklebacks/PLINK/stickleback_all_clean"}

PHENO_DIR=${PHENO_DIR:? "ERROR: You must provide PHENO_DIR via sbatch --export=PHENO_DIR=..."}
QCOVAR=${QCOVAR:? "ERROR: You must provide QCOVAR via sbatch --export=QCOVAR=..."}
COVAR=${COVAR:-"/gpfs01/home/shxyt1/data/sticklebacks/pheno/traits/sex.txt"}

OUTROOT=${OUTROOT:? "ERROR: You must provide OUTROOT via sbatch --export=OUTROOT=..."}

# GRM will always be in OUTROOT
GRM_DIR="${OUTROOT}/GRM"
mkdir -p "${GRM_DIR}"

echo "--------------------------------------"
echo "PHENO_DIR = ${PHENO_DIR}"
echo "QCOVAR = ${QCOVAR}"
echo "COVAR = ${COVAR}"
echo "OUTROOT = ${OUTROOT}"
echo "GRM_DIR = ${GRM_DIR}"
echo "--------------------------------------"

# ============================================================
# 2. 构建 GRM （只做一次）
# ============================================================
if [[ -f "${GRM_DIR}/stickleback_clean.grm.bin" ]]; then
    echo "[GRM] GRM already exists — skip building."
else
    echo "[GRM] Building GRM..."
    ${GCTA_BIN} \
        --bfile "${BFILE}" \
        --make-grm \
        --thread-num 8 \
        --out "${GRM_DIR}/stickleback_clean"
fi

# ============================================================
# 3. 获取所有 phenotype 文件（自动识别 pheno_*.txt）
# ============================================================
PHENO_FILES=($(ls ${PHENO_DIR}/pheno_*.txt | grep -P 'pheno_(As|Mg|Co|Cd|Sr|Pb|K|Fe|Ca)_no_outlier\.txt|pheno_(P|V)_log\.txt' | sort))

TOTAL=${#PHENO_FILES[@]}

echo "--------------------------------------"
echo "Total phenotypes to process: ${TOTAL}"
if [[ ${TOTAL} -gt 0 ]]; then
    echo "First file: ${PHENO_FILES[0]}"
    echo "Last file:  ${PHENO_FILES[-1]}"
fi
echo "--------------------------------------"

if [[ ${SLURM_ARRAY_TASK_ID} -gt ${TOTAL} ]]; then
    echo "Error: array index ${SLURM_ARRAY_TASK_ID} > number of phenotype files (${TOTAL})"
    exit 1
fi
# PHENO_FILES=(
#   ${PHENO_DIR}/pheno_P_log.txt
#   ${PHENO_DIR}/pheno_V_log.txt
# )

TOTAL=${#PHENO_FILES[@]}

echo "Total phenotypes: ${TOTAL}"

if [[ ${SLURM_ARRAY_TASK_ID} -gt ${TOTAL} ]]; then
    echo "Error: array index ${SLURM_ARRAY_TASK_ID} > number of phenotype files (${TOTAL})"
    exit 1
fi

PHENO=${PHENO_FILES[$((SLURM_ARRAY_TASK_ID-1))]}

base=$(basename "${PHENO}")
element=$(echo "${base}" | sed 's/pheno_//; s/.txt//')

echo "Processing phenotype file: ${PHENO}"
echo "Element: ${element}"

OUTDIR="${OUTROOT}/${element}"
mkdir -p "${OUTDIR}"

# ============================================================
# Step 1: REML
# ============================================================
echo "[1] Running REML..."
${GCTA_BIN} \
    --grm "${GRM_DIR}/stickleback_clean" \
    --pheno "${PHENO}" \
    --qcovar "${QCOVAR}" \
    --covar "${COVAR}" \
    --reml \
    --thread-num 8 \
    --out "${OUTDIR}/${element}_REML"

# ============================================================
# Step 2: MLMA GWAS
# ============================================================
echo "[2] Running MLMA GWAS..."
${GCTA_BIN} \
    --bfile "${BFILE}" \
    --grm "${GRM_DIR}/stickleback_clean" \
    --pheno "${PHENO}" \
    --qcovar "${QCOVAR}" \
    --covar "${COVAR}" \
    --mlma \
    --thread-num 8 \
    --out "${OUTDIR}/${element}_MLMA"

echo "Finished ${element} at: $(date)"

# sbatch --export=PHENO_DIR=/gpfs01/home/shxyt1/data/sticklebacks/pheno/traits/log_no_outlier,\
# QCOVAR=/gpfs01/home/shxyt1/data/sticklebacks/pheno/traits/pc5_eda_10kb_hetro_class.txt,\
# OUTROOT=/gpfs01/home/shxyt1/data/sticklebacks/gwas/gwas_5pc_log_EDA10 \
# /gpfs01/home/shxyt1/code/stickleback_qtl_ions/14.3GWAS_log_eda.sh

