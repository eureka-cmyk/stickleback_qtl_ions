#!/bin/bash
#SBATCH --job-name=vcf2gwas_Na
#SBATCH --partition=defq
#SBATCH --mem=100G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/vcf2gwas_Na_%j.log

set -euo pipefail

# 你 HPC 没有 vcf2gwas module，需要 conda（按你老师脚本逻辑）
source $HOME/.bash_profile
conda activate GWAS

VCF="/gpfs01/home/shxyt1/data/sticklebacks/VCF/merge/stickleback_allchr_full.vcf.gz"
PHENO="/gpfs01/home/shxyt1/data/sticklebacks/pheno/vcf2gwas_pheno_Na.tsv"
OUTDIR="/gpfs01/home/shxyt1/data/sticklebacks/analyses/vcf2gwas_Na"

mkdir -p "${OUTDIR}"

# 先不加 GFF：你目前没有注释文件，先跑通统计流程（与GWAS本体无关）
vcf2gwas \
-v "${VCF}" \
-pf "${PHENO}" \
-pheno Na \
-lmm \
-o "${OUTDIR}" \
--memory 99000 \
--threads 20

conda deactivate
