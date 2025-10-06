#!/usr/bin/env python3
from collections import Counter
import os
from pyfaidx import Fasta
import pandas as pd
from Bio.Align import PairwiseAligner


# --- Configuration ---
RAW_DIR = "data/raw"
OUTPUT_DIR = "data/processed"
MARGIN = 2000
GENOME_FASTA = RAW_DIR + "/hg38.fa"
INPUT_BED = RAW_DIR + "/ERVmap.bed"
RMSK_FILE = RAW_DIR + "/rmsk.txt"
OUTPUT_FILE = OUTPUT_DIR + "/paired_ltr_results.txt"

try:
    os.mkdir(OUTPUT_DIR)
except FileExistsError:
    pass

genome = Fasta(GENOME_FASTA)

# Load RepeatMasker
print("Loading RepeatMasker data")
rmsk_cols = [
    "bin", "swScore", "milliDiv", "milliDel", "milliIns",
    "genoName", "genoStart", "genoEnd", "genoLeft",
    "strand", "repName", "repClass", "repFamily",
    "repStart", "repEnd", "repLeft", "id"
]
rmsk = pd.read_csv(RMSK_FILE, sep="\t", header=None, names=rmsk_cols)
rmsk["genoStart"] = rmsk["genoStart"].astype(int)
rmsk["genoEnd"] = rmsk["genoEnd"].astype(int)

print("Filtering RepeatMasker data")
rmsk = rmsk[(rmsk["repClass"] == "LTR") & ((rmsk["repName"].str.contains("LTR")) | (rmsk["repName"].str.contains("-int")))]

aligner = PairwiseAligner()
aligner.mode = 'global'
aligner.match_score = 1
aligner.mismatch_score = 0
aligner.open_gap_score = -5
aligner.extend_gap_score = -0.5

def get_ltr_pair(repNames):
    c = Counter(repNames)
    result = None
    for key in c:
        if '-int' in key:
            continue
        if c[key] > 2:
            return None
        if c[key] == 2:
            if result:
                return None
            else:
                result = key
    if not result:
        return None
    start = repNames.index(result)
    i = start + 1
    has_internal = False
    while repNames[i] != result:
        if '-int' not in repNames[i]:
            return None
        else:
            has_internal = True
        i += 1
    if not has_internal:
        return None
    return result, start, i


def calc_single_point_similarity(alignment):
    diff_counter = 0
    len_counter = 0
    a_gap = False
    b_gap = False
    for i in range(len(alignment[0])):
        a = alignment[0][i]
        b = alignment[1][i]
        if a == b:
            len_counter += 1
            a_gap = False
            b_gap = False
        elif a == '-':
            if not a_gap:
                diff_counter += 1
                len_counter += 1
                a_gap = True
            b_gap = False
        elif b == '-':
            if not b_gap:
                diff_counter += 1
                len_counter += 1
                b_gap = True
            a_gap = False
        else:
            diff_counter += 1
            len_counter += 1
    return (len_counter - diff_counter) / max(len_counter, 1)
  

def process_erv_line(chrom, start, end):
    extended_start = max(0, start - MARGIN)
    extended_end = end + MARGIN

    # LTR-annotations inside range
    subset = rmsk[(rmsk["genoName"] == chrom) &
                  (rmsk["genoStart"] < extended_end) &
                  (rmsk["genoEnd"] > extended_start)]

    if subset.empty:
        return None
    ltr_pair_data = get_ltr_pair(list(subset["repName"]))
    if not ltr_pair_data:
        return None
    ltr_name, first_index, second_index = ltr_pair_data

    # Getting sequences
    start1 = subset.iloc[first_index]["genoStart"]
    end1 = subset.iloc[first_index]["genoEnd"]
    start2 = subset.iloc[second_index]["genoStart"]
    end2 = subset.iloc[second_index]["genoEnd"]
    first_seq = genome[chrom][start1:end1].seq
    second_seq = genome[chrom][start2:end2].seq
    alignments = aligner.align(first_seq, second_seq)
    best_alignment = alignments[0]
    if not first_seq or not second_seq:
        return None
    similarity = calc_single_point_similarity(best_alignment)
    best_hit = (
        chrom,
        f"{chrom}:{start1}-{end1}",
        f"{chrom}:{start2}-{end2}",
        end1 - start1,
        end2 - start2,
        similarity,
        ltr_name
    )
    return best_hit

# --- Main loop ---
print("Processing " + INPUT_BED)
with open(INPUT_BED, 'r') as bed, open(OUTPUT_FILE, 'w') as out_f:
    for line in bed:
        parts = line.strip().split('\t')
        chrom, start, end = parts[0], int(parts[1]), int(parts[2])
        chrom = 'chr' + chrom
        hit = process_erv_line(chrom, start, end)
        if hit:
            out_f.write(f"{hit[0]} {hit[1]} {hit[2]} {hit[3]} {hit[4]} {hit[5]:.3f} {hit[6]}\n")

print(f"Done. Results are in {OUTPUT_FILE}")

