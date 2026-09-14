#!/bin/bash
#SBATCH --job-name=make_pheno_files
#SBATCH --partition=defq
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs01/home/shxyt1/slurm_outputs/make_pheno_files_%j.log

module load python-uoneasy/3.12.3-GCCcore-13.3.0

PHENO_CSV="/gpfs01/home/shxyt1/data/sticklebacks/pheno/qtl_all_phenotype_data.csv"
OUTDIR="/gpfs01/home/shxyt1/data/sticklebacks/pheno"

mkdir -p "$OUTDIR"

python << 'EOF'
import pandas as pd
import os

pheno_csv = "/gpfs01/home/shxyt1/data/sticklebacks/pheno/qtl_all_phenotype_data.csv"
outdir = "/gpfs01/home/shxyt1/data/sticklebacks/pheno"

df = pd.read_csv(pheno_csv)

# rename ID column to FID/IID for PLINK format
if "fish_id" in df.columns:
    df.rename(columns={"fish_id": "IID"}, inplace=True)
df["FID"] = df["IID"]

# columns not to export as phenotype
exclude = {
    "IID", "FID",
    "date_fish_caught",
    "date_fish_processed",
    "sex",
    "notes",
    "icpms_run",
    "schisto"
}

traits = [c for c in df.columns if c not in exclude]

print(f"Found {len(traits)} phenotype traits:")
for t in traits:
    print(" -", t)

# export each phenotype to its own text file
for trait in traits:
    outpath = os.path.join(outdir, f"pheno_{trait}.txt")
    sub = df[["FID", "IID", trait]].copy()
    sub.to_csv(outpath, sep="\t", index=False, header=True)

print("\nAll phenotype files written to:", outdir)
EOF
