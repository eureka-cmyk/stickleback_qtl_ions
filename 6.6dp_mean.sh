#!/bin/bash
#SBATCH --job-name=meanDP
#SBATCH --partition=defq
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --array=1-3

DP_DIR="/share/MacColl_shared/QTL/vcf/GCF_016920845.1_GAculeatus_UGA_version5_genomic/sex_coverage_DP"
VCF_DIR="/share/MacColl_shared/QTL/vcf/GCF_016920845.1_GAculeatus_UGA_version5_genomic/ploidy_aware_HWEPops_MQ10_BQ20"
OUT_DIR="$DP_DIR/mean2"
mkdir -p $OUT_DIR

CHROMS=("NC_053212.1" "NC_053230.1" "NC_053233.1")
CHR=${CHROMS[$SLURM_ARRAY_TASK_ID-1]}

DP_FILE="$DP_DIR/${CHR}_DP.txt"
VCF_FILE="$VCF_DIR/stickleback_${CHR}_sorted.bcf"
OUT="$OUT_DIR/${CHR}_mean.txt"

echo "Processing $CHR"

############################################
# 1. Get sample IDs in correct order
############################################
bcftools query -l "$VCF_FILE" > "$OUT.samples"

############################################
# 2. Compute mean DP for each sample column
############################################
awk '
NR==1 {
    # number of samples = NF - 2 (CHR,POS)
    for(i=3;i<=NF;i++){
        sum[i]=0
        count[i]=0
    }
}
NR>1 {
    for(i=3;i<=NF;i++){
        sum[i] += $i
        count[i]++
    }
}
END{
    for(i=3;i<=NF;i++){
        if(count[i] > 0)
            print sum[i] / count[i]
        else
            print 0
    }
}' "$DP_FILE" > "$OUT.values"

############################################
# 3. Join sampleID + meanDP
############################################
paste "$OUT.samples" "$OUT.values" > "$OUT"

rm "$OUT.samples" "$OUT.values"

echo "Done: $OUT"
