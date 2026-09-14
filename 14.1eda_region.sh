VCF="/gpfs01/home/shxyt1/data/sticklebacks/VCF/merge/stickleback_allchr_full.vcf.gz"
REGION2="NC_053215.1:12683708-12893934"
REGION_NEW="NC_053215.1:12783500-12800300"
OUT2="eda_NC_053215.1_12683708_12893934.vcf.gz"
OUT_PREFIX="eda_subset_16kb"
OUT_VCF="${OUT_PREFIX}.vcf.gz"
OUT_BISNP="${OUT_PREFIX}.bisnp.vcf.gz"

module load bcftools-uoneasy/1.19-GCC-13.2.0

bcftools view -r "$REGION2" -Oz -o "$OUT2" "$VCF"
bcftools index -t "$OUT2"
bcftools view -H "$OUT2" | wc -l
#209392

# 2. 提取特定区域
bcftools view -r "$REGION_NEW" -Oz -o "$OUT_VCF" "$VCF"
bcftools index -t "$OUT_VCF"

# 3. 过滤出二等位基因 SNP (Biallelic SNPs)
bcftools view \
  -i 'TYPE="snp" && N_ALT=1' \
  -Oz \
  -o "$OUT_BISNP" \
  "$OUT_VCF"


bcftools view \
>   -i 'TYPE="snp" && N_ALT=1' \
>   -Oz \
>   -o eda_NC_053215.1_12683708_12893934.bisnp.vcf.gz \
>   eda_NC_053215.1_12683708_12893934.vcf.gz

bcftools index -t eda_NC_053215.1_12683708_12893934.bisnp.vcf.gz

bcftools view -H eda_NC_053215.1_12683708_12893934.bisnp.vcf.gz | wc -l
#5227

# 4. 转换为 PLINK 格式并进行初步质量控制
module load plink-uoneasy/1.9b_6.21-x86_64

plink \
  --vcf "$OUT_BISNP" \
  --allow-extra-chr \
  --double-id \
  --geno 0.2 \
  --mind 0.1 \
  --make-bed \
  --out "${OUT_PREFIX}.filt"

# plink \
#   --vcf "$OUT_BISNP" \
#   --allow-extra-chr \
#   --double-id \
#   # 几乎关掉所有过滤
#   --geno 0.99 \ 
#   --mind 0.99 \
#   --make-bed \
#   --out "${OUT_PREFIX}.raw_comparison"

plink \
  --vcf eda_NC_053215.1_12683708_12893934.bisnp.vcf.gz \
  --allow-extra-chr \
  --double-id \
  --make-bed \
  --out eda_edaRegion

# 5. 计算杂合度 (Heterozygosity)
plink \
  --bfile "${OUT_PREFIX}.filt" \
  --allow-extra-chr \
  --het \
  --out "${OUT_PREFIX}_het"



# 直接用 filt 文件跑 PCA，跳过 prune.in 提取
# 使用 10kb 区域内未经过 LD 修剪的所有 774 个位点跑 PCA
# plink \
#   --bfile /gpfs01/home/shxyt1/data/sticklebacks/VCF/eda/eda10kb.filt \
#   --allow-extra-chr \
#   --pca 10 \
#   --out /gpfs01/home/shxyt1/data/sticklebacks/VCF/eda/eda10kb_no_prune.pca


ls -lh eda_edaRegion.bed eda_edaRegion.bim eda_edaRegion.fam eda_edaRegion.log

# 6.PCA
# 跳过 extract 步骤，直接对过滤后的 412 个变异位点跑 PCA
plink \
  --bfile eda_subset_16kb.filt \
  --allow-extra-chr \
  --pca 10 \
  --out eda_subset_16kb.pca


# plink \
#   --bfile eda10kb \
#   --allow-extra-chr \
#   --geno 0.2 \
#   --mind 0.1 \
#   --make-bed \
#   --out eda10kb.filt

plink \
  --bfile eda10kb \
  --allow-extra-chr \
  --geno 0.01 \
  --mind 0.1 \
  --make-bed \
  --out eda10kb_geno0.01.filt

plink \
  --bfile eda_edaRegion.filt \
  --allow-extra-chr \
  --indep-pairwise 50 5 0.2 \
  --out eda_edaRegion.prune

ls -lh eda_edaRegion.prune.prune.in eda_edaRegion.prune.prune.out

plink \
  --bfile eda_edaRegion.filt \
  --allow-extra-chr \
  --extract eda_edaRegion.prune.prune.in \
  --pca 10 \
  --out eda_edaRegion.pca

module load plink-uoneasy/1.9b_6.21-x86_64

#hetro
plink \
  --bfile /gpfs01/home/shxyt1/data/sticklebacks/VCF/eda/eda_edaRegion.filt \
  --allow-extra-chr \
  --het \
  --out /gpfs01/home/shxyt1/data/sticklebacks/VCF/eda/eda25kb

#mds
plink --bfile eda10kb.filt \
>       --extract eda10kb.prune.prune.in \
>       --cluster \
>       --mds-plot 4 \
>       --allow-extra-chr \
>       --out eda10kb_mds

#missing

plink --bfile /gpfs01/home/shxyt1/data/sticklebacks/VCF/eda/eda10kb \
      --allow-extra-chr \
      --missing \
      --out sample_missing_10kb



plink --bfile /gpfs01/home/shxyt1/data/sticklebacks/VCF/eda/eda10kb \
      --allow-extra-chr \
      --recode A \
      --out eda10kb_matrix

#name snp
plink \
  --bfile eda10kb_geno0.01.filt \
  --allow-extra-chr \
  --set-missing-var-ids @:#\$1,\$2 \
  --make-bed \
  --out eda10kb_geno0.01.filt.named

plink \
  --bfile eda10kb_geno0.01.filt.named \
  --allow-extra-chr \
  --recode A \
  --out eda10kb_geno0.01.filt.named

#missing
plink \
  --bfile eda10kb_geno0.01.filt.named \
  --allow-extra-chr \
  --missing \
  --out eda10kb_geno0.01.filt.named

#pca
plink \
  --bfile eda10kb_geno0.01.filt.named \
  --allow-extra-chr \
  --pca 10 \
  --out eda10kb_geno0.01.filt.named

#MDS
plink \
  --bfile eda10kb_geno0.01.filt.named \
  --allow-extra-chr \
  --cluster \
  --mds-plot 4 \
  --out eda10kb_geno0.01.filt.named

#sample num
wc -l eda10kb_geno0.1.filt.named.fam

#snp num
wc -l eda10kb_geno0.01.filt.bim