#!/usr/bin/env bash
set -euo pipefail

# Work from repo root
BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")"/.. && pwd)"
cd "$BASE_DIR"

RAW_DIR="./data/raw/maf"
OUT_DIR="./data/processed/maf_filtered"
SCRIPT="./src/filter_maf.py"

mkdir -p "$OUT_DIR"

echo "[INFO] Filtering MAF files..."
CHROMS=$(ls "$RAW_DIR" | sed -E 's/\.maf\.gz//' | tr '\n' ' ')
echo "[INFO] Found $(echo "$CHROMS" | wc -w) chromosomes."

# Run in parallel with visible progress
parallel -j 10 "python3 $SCRIPT {1} && echo '[DONE] {1}'" ::: $CHROMS

echo "[INFO] All filtered MAFs saved to $OUT_DIR"

