#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot PCA (step4_PCA.eigenvec) and Sample Missingness (all_samples_missing.imiss)
Yiwen Tao – GWAS QC pipeline
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

##############################
# 1. FILE PATHS
##############################

# Automatically detect working directory (plink_check)
workdir = "/gpfs01/home/shxyt1/data/sticklebacks/plink_check"

pca_file = f"{workdir}/step4_PCA.eigenvec"
missing_file = f"{workdir}/all_samples_missing.imiss"

out_pca_png = f"{workdir}/pca_plot.png"
out_missing_png = f"{workdir}/missingness_plot.png"


##############################
# 2. PCA PLOT
##############################

print("[Plot] PCA scatter plot...")

# eigenvec columns: FID IID PC1 PC2 ...
pca = pd.read_csv(pca_file, delim_whitespace=True, header=None)
pca.columns = ["FID", "IID"] + [f"PC{i}" for i in range(1, pca.shape[1]-1)]

plt.figure(figsize=(8, 6))
sns.scatterplot(
    x="PC1",
    y="PC2",
    data=pca,
    s=25,
    edgecolor="black",
    linewidth=0.3
)

plt.title("PCA: PC1 vs PC2", fontsize=14)
plt.xlabel("PC1")
plt.ylabel("PC2")

plt.tight_layout()
plt.savefig(out_pca_png, dpi=300)
plt.close()

print(f"[Saved] PCA plot → {out_pca_png}")


##############################
# 3. MISSINGNESS PLOT
##############################

print("[Plot] Sample missingness...")

miss = pd.read_csv(missing_file, delim_whitespace=True)
# imiss columns: FID IID MISS_PHENO N_MISS N_GENO F_MISS

plt.figure(figsize=(10, 4))
sns.barplot(
    x="IID",
    y="F_MISS",
    data=miss,
    color="steelblue"
)
plt.xticks([], [])  # hide sample labels
plt.ylabel("Missing rate (F_MISS)")
plt.title("Per-sample missing genotype rate")

plt.axhline(0.10, color="red", linestyle="--", label="10% cutoff")
plt.legend()

plt.tight_layout()
plt.savefig(out_missing_png, dpi=300)
plt.close()

print(f"[Saved] Missingness plot → {out_missing_png}")

print("All plots generated successfully.")
