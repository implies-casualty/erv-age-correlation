#!/bin/bash
# Download hg38 genome, ERVmap, RepeatMasker rmsk, and multiz30way MAF alignments

set -euo pipefail

# Base dirs
RAW_DIR="data/raw"
MAF_DIR="${RAW_DIR}/maf"

mkdir -p "${RAW_DIR}" "${MAF_DIR}"

#############################################
# 1. Human genome (hg38)
#############################################
echo "Downloading hg38 reference genome..."
wget -c "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz" -O "${RAW_DIR}/hg38.fa.gz"
gunzip -kf "${RAW_DIR}/hg38.fa.gz"   # keep compressed & uncompressed copies

#############################################
# 2. RepeatMasker (rmsk)
#############################################
echo "Downloading RepeatMasker annotations (rmsk.txt.gz)..."
if command -v rsync >/dev/null 2>&1; then
    rsync -avz rsync://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz "${RAW_DIR}/"
else
    echo "rsync not found — falling back to wget"
    wget -c "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz" -O "${RAW_DIR}/rmsk.txt.gz"
fi
gunzip -kf "${RAW_DIR}/rmsk.txt.gz"

#############################################
# 3. ERVmap annotation
#############################################
echo "Downloading ERVmap annotations..."
wget -c "https://raw.githubusercontent.com/mtokuyama/ERVmap/master/ref/ERVmap.bed" -O "${RAW_DIR}/ERVmap.bed"

#############################################
# 4. Multiz30way MAF alignments (per chromosome)
#############################################
echo "Downloading multiz30way MAF alignments (chr1–22, X, Y)..."
BASE_URL="http://hgdownload.soe.ucsc.edu/goldenPath/hg38/multiz30way/maf"

# chr1–22
for i in {1..22}; do
    wget -c "${BASE_URL}/chr${i}.maf.gz" -O "${MAF_DIR}/chr${i}.maf.gz"
done

# chrX, chrY
for c in X Y; do
    wget -c "${BASE_URL}/chr${c}.maf.gz" -O "${MAF_DIR}/chr${c}.maf.gz"
done

echo "✅ All downloads complete. Files saved in ${RAW_DIR}/"

