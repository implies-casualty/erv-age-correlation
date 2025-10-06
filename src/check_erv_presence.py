#!/usr/bin/env python3
"""
Determine which primate genomes share each human ERV locus by analyzing
filtered multiz30way MAF alignments. For every human LTR pair, the script
measures how much of its region aligns (non-gap bases) in other species.

Input:
    data/processed/paired_ltr_results.txt   # list of LTR pairs and coordinates
    data/processed/maf_filtered/*.ERV.maf.gz
Output:
    data/processed/ltr_presence.csv
"""

from collections import defaultdict
import gzip
from intervaltree import IntervalTree
import os
import pandas as pd

from Bio import AlignIO

# ---------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------
PROCESSED_DIR = "./data/processed"
bed_file = PROCESSED_DIR + "/paired_ltr_results.txt"
maf_dir = PROCESSED_DIR + "/maf_filtered"
out_file = PROCESSED_DIR + "/ltr_presence.csv"

# ---------------------------------------------------------------------
# Step 1: Load LTR coordinates
# ---------------------------------------------------------------------
print("[INFO] Loading LTR coordinates...")
ltr_regions = defaultdict(IntervalTree)

with open(bed_file) as f:
    for line in f:
        if not line.strip():
            continue
        chrom, loc1, loc2, _ = line.strip().split(" ", 3)
        for loc in (loc1, loc2):
            chrom_name, start_end = loc.split(":", 1)
            chrom = chrom_name.replace("chr", "")
            start, end = map(lambda x: int(float(x)), start_end.split("-"))
            ltr_regions[chrom].addi(start, end, loc)

print(f"[INFO] Loaded {sum(len(t) for t in ltr_regions.values())} LTR regions.")

# ---------------------------------------------------------------------
# Step 2: Iterate through MAFs and calculate coverage
# ---------------------------------------------------------------------
agg_results = defaultdict(int)
ltr_lengths = {}

for chrom in ltr_regions:
    maf_file = maf_dir + f"/chr{chrom}.ERV.maf.gz"
    if not os.path.exists(maf_file):
        print(f"[WARN] Missing {maf_file}")
        continue

    print(f"[INFO] Processing {maf_file} ...")
    with gzip.open(maf_file, "rt") as handle:
        for alignment in AlignIO.parse(handle, "maf"):
            # Locate human component
            human_seq = next((r for r in alignment if r.id.startswith("hg38.")), None)
            if human_seq is None:
                continue

            human_start = human_seq.annotations["start"]
            human_end = human_start + human_seq.annotations["size"]

            # Check overlaps with LTR intervals
            for ov in ltr_regions[chrom].overlap(human_start, human_end):
                ltr_name = ov.data
                ltr_len = ov.end - ov.begin
                ltr_lengths[(chrom, ov.begin, ov.end, ltr_name)] = ltr_len

                # Map alignment columns to genome coordinates
                aln_to_genome = []
                gpos = human_start
                for base in str(human_seq.seq):
                    if base != "-":
                        aln_to_genome.append(gpos)
                        gpos += 1
                    else:
                        aln_to_genome.append(None)

                # LTR indices in this alignment block
                ltr_indices = [
                    i for i, g in enumerate(aln_to_genome)
                    if g is not None and ov.begin <= g < ov.end
                ]
                if not ltr_indices:
                    continue

                # Measure aligned coverage in other species
                for rec in alignment:
                    species = rec.id.split(".")[0]
                    if species == "hg38":
                        continue
                    subseq = "".join(rec.seq[i] for i in ltr_indices)
                    aligned_bases = sum(1 for b in subseq if b != "-")
                    key = (chrom, ov.begin, ov.end, ltr_name, species)
                    agg_results[key] += aligned_bases

# ---------------------------------------------------------------------
# Step 3: Build summary table
# ---------------------------------------------------------------------
print("[INFO] Building output table...")
records = []
for key, aligned_bases in agg_results.items():
    chrom, start, end, ltr_name, species = key
    ltr_len = ltr_lengths.get((chrom, start, end, ltr_name), 0)
    coverage = aligned_bases / ltr_len if ltr_len > 0 else 0
    records.append({
        "chrom": chrom,
        "start": start,
        "end": end,
        "ltr": ltr_name,
        "species": species,
        "aligned_bases": aligned_bases,
        "ltr_len": ltr_len,
        "coverage": round(coverage, 3),
    })

df = pd.DataFrame(records)
df.to_csv(out_file, index=False)
print(f"[INFO] Results written to {out_file}")

