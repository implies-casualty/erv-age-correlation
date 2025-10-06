#!/usr/bin/env python3
import sys
import gzip
from bx.align import maf
from bx.intervals import IntervalTree
from collections import defaultdict

MARGINS = 2000

RAW_DIR = "data/raw"
OUTPUT_DIR = "data/processed"

bed_file = RAW_DIR + "/ERVmap.bed"
maf_dir = RAW_DIR + "/maf"
out_dir = OUTPUT_DIR + "/maf_filtered"

# --- Load all BED intervals grouped by chromosome ---
intervals_by_chr = defaultdict(IntervalTree)
with open(bed_file) as bed:
    for line in bed:
        if line.startswith("#") or not line.strip():
            continue
        chrom, start, end, *_ = line.strip().split()
        # strip "chr" prefix if present
        chrom = chrom.replace("chr", "")
        intervals_by_chr[chrom].add(int(start) - MARGINS, int(end) + MARGINS)


def process_chr(chrom):
    in_maf = f"{maf_dir}/{chrom}.maf.gz"
    out_maf = f"{out_dir}/{chrom}.ERV.maf.gz"
    tree = intervals_by_chr.get(chrom.replace("chr", ""), None)

    if tree is None:
        return

    with gzip.open(in_maf, "rt") as inf, gzip.open(out_maf, "wt") as outf:
        reader = maf.Reader(inf)
        writer = maf.Writer(outf)

        for m in reader:
            comp = m.components[0]
            block_start = comp.start
            block_end = comp.start + comp.size
            if tree.find(block_start, block_end):
                writer.write(m)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: genome_filter_maf.py chrN")
        sys.exit(1)
    process_chr(sys.argv[1])

