#!/bin/bash
#SBATCH --job-name=install_gcta
#SBATCH --partition=defq
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=gcta_install_%j.log

echo "=== Starting GCTA installation ==="

# Installation directory
INSTALL_DIR=/gpfs01/home/shxyt1/app/gcta


# Choose GCTA version
GCTA_VERSION="1.95"
GCTA_ZIP="gcta_${GCTA_VERSION}.zip"
DOWNLOAD_URL="https://yanglab.westlake.edu.cn/software/gcta/bin/gcta-1.95.0-linux-kernel-3-x86_64.zip"

echo "Downloading GCTA from:"
echo $DOWNLOAD_URL

# Download
cd $INSTALL_DIR
wget -O $GCTA_ZIP $DOWNLOAD_URL

# unzip gcta_1.95.zip
unzip $GCTA_ZIP
rm $GCTA_ZIP

# Make executable
chmod +x gcta64

echo "=== GCTA installed successfully at: $INSTALL_DIR/gcta64 ==="

echo "Add the following line to your ~/.bashrc:"
echo "    export PATH=\$PATH:${INSTALL_DIR}"

export PATH=\$PATH:$/gpfs01/home/shxyt1/app/gcta/gcta-1.95.0-linux-kernel-3-x86_64
