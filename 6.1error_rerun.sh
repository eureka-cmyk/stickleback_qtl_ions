#!/bin/bash
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=40G
#SBATCH --time=12:00:00
#SBATCH --job-name=post-snp
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/sort-snp-%A.out

# === 路径 ===
genome_name=GCF_016920845.1_GAculeatus_UGA_version5_genomic
master_output=/share/MacColl_shared/QTL/vcf/${genome_name}/ploidy_aware_HWEPops_MQ10_BQ20
BED=/share/MacColl_shared/QTL/genomes/test.bed

# === 关键修复：指定大临时目录（而不是 /tmp）===
TMPDIR=/share/MacColl_shared/QTL/tmp_sort
mkdir -p ${TMPDIR}

# 读取染色体名称
chroms=$(awk '{print $1}' $BED)

for chr in $chroms; do
    echo "============================================"
    echo "Processing chromosome: $chr"
    echo "============================================"

    VCF_RAW=${master_output}/stickleback_${chr}.bcf
    VCF_SORTED=${master_output}/stickleback_${chr}_sorted.bcf
    VCF_SORTED_INDEX=${VCF_SORTED}.csi

    if [ ! -s "$VCF_RAW" ]; then
        echo "[SKIP] Raw BCF not found or empty: $VCF_RAW"
        continue
    fi

    # 为该染色体创建独立临时目录
    LOCAL_TMP=${TMPDIR}/sort_${chr}
    mkdir -p ${LOCAL_TMP}
    echo "[INFO] Sorting with TMPDIR=${LOCAL_TMP}"

    bcftools sort \
        -T ${LOCAL_TMP} \
        -O b \
        -o $VCF_SORTED \
        $VCF_RAW

    echo "[INFO] Indexing..."
    bcftools index $VCF_SORTED

    if [ -f "$VCF_SORTED_INDEX" ]; then
        echo "[INFO] Removing unsorted BCF"
        rm $VCF_RAW
    fi

    echo "[INFO] Counting variants..."
    bcftools query -l $VCF_SORTED | wc -l | xargs echo "Samples:"
    bcftools view -H $VCF_SORTED | wc -l | xargs echo "Variants:"
    bcftools view -H -v snps $VCF_SORTED | wc -l | xargs echo "SNPs:"
    bcftools view -H -v indels $VCF_SORTED | wc -l | xargs echo "Indels:"

    echo "[DONE] Chromosome $chr finished."
    echo
done

# # === 染色体名称 ===
# #chr=$(awk "NR==$SLURM_ARRAY_TASK_ID" /share/MacColl_shared/QTL/genomes/test.bed | awk '{print $1}')

# VCF_RAW=${master_output}/stickleback_${chr}.bcf
# VCF_SORTED=${master_output}/stickleback_${chr}_sorted.bcf
# VCF_SORTED_INDEX=${master_output}/stickleback_${chr}_sorted.bcf.csi


# echo "Processing chromosome: $chr"
# echo " RAW:    $VCF_RAW"
# echo " SORTED: $VCF_SORTED"

# # === raw bcf 是否存在 ===
# if [ ! -s "$VCF_RAW" ]; then
#     echo "[ERROR] RAW BCF does not exist or is empty. Skip."
#     exit 1
# fi

# # === 排序 ===
# echo "[INFO] Sorting..."

# # 设置一个大 tmpdir
# TMPDIR=/share/MacColl_shared/QTL/tmp_sort
# mkdir -p ${TMPDIR}/sort_${chr}

# bcftools sort \
#     -T ${TMPDIR}/sort_${chr} \
#     -O b \
#     -o $VCF_SORTED \
#     $VCF_RAW

# #bcftools sort -O b -o $VCF_SORTED $VCF_RAW

# # === 建 index ===
# echo "[INFO] Indexing..."
# bcftools index $VCF_SORTED

# # === 删除原始未排序的 bcf ===
# if [ -f "$VCF_SORTED_INDEX" ]; then
#     echo "[INFO] Removing unsorted BCF: $VCF_RAW"
#     rm $VCF_RAW
# fi

# echo "[INFO] Sorted BCF created: $VCF_SORTED"

# # === 输出统计 ===
# echo "These are the individuals in the VCF file:"
# bcftools query -l $VCF_SORTED

# echo "This is the number of variants in the file:"
# bcftools view -H $VCF_SORTED | wc -l

# echo "This is the number of SNPs in the file:"
# bcftools view -H -v snps $VCF_SORTED | wc -l

# echo "This is the number of indels in the file:"
# bcftools view -H -v indels $VCF_SORTED | wc -l

# echo "[DONE] Chromosome $chr processing complete."
