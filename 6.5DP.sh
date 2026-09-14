#!/bin/bash
#SBATCH --job-name=extractDP
#SBATCH --output=extractDP_%j.out
#SBATCH --error=extractDP_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --partition=defq

module load bcftools-uoneasy/1.19-GCC-13.2.0

# 输入 BCF 目录
BCF_DIR="/share/MacColl_shared/QTL/vcf/GCF_016920845.1_GAculeatus_UGA_version5_genomic/ploidy_aware_HWEPops_MQ10_BQ20"

# 输出目录
OUT_DIR="/share/MacColl_shared/QTL/vcf/GCF_016920845.1_GAculeatus_UGA_version5_genomic/sex_coverage_DP"
mkdir -p "$OUT_DIR"

# 要提取 DP 的三个染色体
CHROMS=("NC_053212.1" "NC_053230.1" "NC_053233.1")

echo "=== Extracting DP for 3 chromosomes ==="
echo "Output directory: $OUT_DIR"
echo ""

for chr in "${CHROMS[@]}"; do
    BCF_FILE="$BCF_DIR/stickleback_${chr}_sorted.bcf"
    OUT_FILE="$OUT_DIR/${chr}_DP.txt"

    echo "Processing $chr ..."
    echo "BCF:  $BCF_FILE"
    echo "OUT:  $OUT_FILE"

    if [[ ! -f "$BCF_FILE" ]]; then
        echo "ERROR: BCF not found: $BCF_FILE"
        exit 1
    fi

    # %DP for each sample → each column is a sample's depth
    bcftools query -f '%CHROM %POS[%DP ]\n' "$BCF_FILE" > "$OUT_FILE"

    echo "Finished $chr"
    echo ""
done

echo "=== ALL DONE ==="
echo "DP files saved in: $OUT_DIR"
