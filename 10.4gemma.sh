#!/bin/bash
#SBATCH --job-name=Na_RINT_BSLMM_simple
#SBATCH --partition=defq
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/Na_RINT_BSLMM_simple_%j.log


# 1. 设置路径
WORKDIR="/gpfs01/home/shxyt1/data/sticklebacks/gemma/Na_RINT_simple"

GENO_PREFIX="/gpfs01/home/shxyt1/data/sticklebacks/PLINK/stickleback_all_clean"
PHENO_FILE="/gpfs01/home/shxyt1/data/sticklebacks/pheno/traits/pheno_Na_RINT.txt"

mkdir -p "${WORKDIR}"
cd "${WORKDIR}"

echo "Start GEMMA BSLMM (simple) for Na_RINT"
date
hostname

# 2. 软链接 PLINK 文件（避免重复拷贝大文件）
ln -sf ${GENO_PREFIX}.bed ./geno.bed
ln -sf ${GENO_PREFIX}.bim ./geno.bim
ln -sf ${GENO_PREFIX}.fam ./geno.fam

# 3. 复制 phenotype 到当前目录（只是方便）
cp "${PHENO_FILE}" ./pheno.txt

# 4. 跑最基础的 BSLMM（不加 covar）
gemma \
  -bfile geno \
  -bslmm 1 \
  -n 1 \
  -maf 0.01 \
  -o Na_RINT_BSLMM_simple \
  -p pheno.txt

echo "Done. Check:"
echo "  ${WORKDIR}/Na_RINT_BSLMM_simple.log.txt"
echo "  ${WORKDIR}/Na_RINT_BSLMM_simple.hyp.txt"
