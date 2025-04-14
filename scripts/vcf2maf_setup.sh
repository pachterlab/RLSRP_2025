#!/bin/bash

set -e  # Exit on any error
set -u  # Treat unset variables as errors

# ================== USER CONFIG ==================
RELEASE="93"  # This entire script assumes homo_sapiens
BUILD="GRCh37"  # or "GRCh38";
MINICONDA_PATH="$HOME/miniconda3"
ENV_NAME="vep${RELEASE}"
VEP_DATA="$HOME/.vep"
OPT_DIR="$HOME"/opt
# =================================================

# Derived
CACHE_FILE="homo_sapiens_vep_${RELEASE}_${BUILD}.tar.gz"
CACHE_URL="ftp://ftp.ensembl.org/pub/release-${RELEASE}/variation/indexed_vep_cache/${CACHE_FILE}"
VEP_VERSION="${RELEASE}"

echo "=== Downloading vcf2maf ==="
mkdir -p "$OPT_DIR"
cd "$OPT_DIR"
VCF2MAF_DIR=$(find "$OPT_DIR" -maxdepth 1 -type d -name "mskcc-vcf2maf-*")
if [ -z "$VCF2MAF_DIR" ]; then
    echo "Downloading and extracting vcf2maf..."
    VCF2MAF_URL=$(curl -sL https://api.github.com/repos/mskcc/vcf2maf/releases | grep -m1 tarball_url | cut -d\" -f4)
    curl -L -o mskcc-vcf2maf.tar.gz "$VCF2MAF_URL"
    tar -zxf mskcc-vcf2maf.tar.gz
    rm mskcc-vcf2maf.tar.gz
    VCF2MAF_DIR=$(find "$OPT_DIR" -maxdepth 1 -type d -name "mskcc-vcf2maf-*" | head -n 1)
else
    echo "vcf2maf already exists: $VCF2MAF_DIR"
fi
cd "$VCF2MAF_DIR"

# Create conda env if not already created
echo "=== Setting up conda environment ==="
source "${MINICONDA_PATH}/etc/profile.d/conda.sh"
if ! conda env list | grep -q "^${ENV_NAME}"; then
    echo "Creating conda env ${ENV_NAME}..."
    conda create -y -n "$ENV_NAME"
    conda activate "$ENV_NAME"
    conda install -y -c conda-forge -c bioconda -c defaults \
        ensembl-vep="$VEP_VERSION" htslib bcftools samtools \
        ucsc-liftover perl-list-moreutils
else
    echo "Conda env ${ENV_NAME} already exists. Activating..."
    conda activate "$ENV_NAME"
fi

# Download VEP cache if not already present
echo "=== Checking for VEP cache ==="
mkdir -p "$VEP_DATA"
if [ ! -d "$VEP_DATA/homo_sapiens/${RELEASE}_${BUILD}" ]; then
    echo "Downloading VEP cache for ${BUILD} release ${RELEASE}..."
    wget "$CACHE_URL"
    tar -xzf "$CACHE_FILE" -C "$VEP_DATA"
else
    echo "VEP cache already exists."
fi

echo "Finished Script"