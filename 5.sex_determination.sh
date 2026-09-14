#!/bin/bash
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --job-name=sexcov
#SBATCH --array=1-352
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/sex-cov-%j.out

module load samtools-uoneasy/1.18-GCC-12.3.0

# 目录设置
BAMDIR="/share/MacColl_shared/QTL/bams/GCF_016920845.1_GAculeatus_UGA_version5_genomic/clean_bams"
OUTDIR="/share/MacColl_shared/QTL/bams/GCF_016920845.1_GAculeatus_UGA_version5_genomic/Sex_coverage"
mkdir -p $OUTDIR

# 获取 bam 列表
BAMLIST=($(ls $BAMDIR/*.bam))
BAM=${BAMLIST[$SLURM_ARRAY_TASK_ID-1]}

# 个体 ID 从文件名提取
IND=$(basename $BAM .bam)

# 输出文件
OUTFILE="$OUTDIR/${IND}_coverage.txt"

echo "Processing $IND"

# 染色体
CHR_AUTO="NC_053212.1"   # 常染色体 
CHR_X="NC_053230.1"      # X
CHR_Y="NC_053233.1"      # Y

# 计算覆盖度
# samtools depth -a \
#     -r $CHR_AUTO \
#     -r $CHR_X \
#     -r $CHR_Y \
#     $BAM \
#     | awk -v IND=$IND '{cov[$1]+=$3; n[$1]++}
#         END{
#             for (c in cov) print IND, c, cov[c]/n[c]
#         }' > $OUTFILE


# 一个函数：给定染色体名，输出平均 depth
mean_depth () {
    samtools depth -a -r "$1" "$BAM" \
    | awk '{s+=$3; n++} END{if(n>0) print s/n; else print 0}'
}

AUTO_MEAN=$(mean_depth $CHR_AUTO)
X_MEAN=$(mean_depth $CHR_X)
Y_MEAN=$(mean_depth $CHR_Y)

# 避免除以 0
if [ "$AUTO_MEAN" = "0" ]; then
    X_RATIO=0
    Y_RATIO=0
else
    X_RATIO=$(awk -v x="$X_MEAN" -v a="$AUTO_MEAN" 'BEGIN{print x/a}')
    Y_RATIO=$(awk -v y="$Y_MEAN" -v a="$AUTO_MEAN" 'BEGIN{print y/a}')
fi

{
echo -e "$IND\t$CHR_AUTO\t$AUTO_MEAN"
echo -e "$IND\t$CHR_X\t$X_MEAN"
echo -e "$IND\t$CHR_Y\t$Y_MEAN"
echo -e "$IND\tX_ratio\t$X_RATIO"
echo -e "$IND\tY_ratio\t$Y_RATIO"
} > "$OUTFILE"