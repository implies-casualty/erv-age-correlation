#!/usr/bin/env python3
"""
summarize_erv_similarity.py

Combine paired LTR similarity data and cross-species ERV presence data to
summarize, for each ERV locus:

  • how many paired LTRs were found nearby
  • the most distant primate species sharing this ERV
  • estimated evolutionary distance (rank and MYA)
  • average similarity of the paired LTRs

Outputs: 'erv_ltr_presence_summary.csv'
"""

import pandas as pd

# ------------------- CONFIGURATION -------------------
ERV_BED = "data/raw/ERVmap.bed"           # Path to ERVmap.bed
PAIRED_LTR = "data/processed/paired_ltr_results.txt"       # Path to paired LTR similarity data
LTR_PRESENCE = "data/processed/ltr_presence.csv"           # Path to species presence/coverage data

OUTPUT = "data/processed/erv_ltr_presence_summary.csv"

MARGIN = 2000                               # Window around ERV to match LTRs (bp)
COVERAGE_THRESHOLD = 0.1                    # LTR considered present if coverage > this

# ------------------- SPECIES LCA DATA -------------------
# Each species has: common_name, estimated LCA in MYA, and a "rank" where
# smaller = closer to human, larger = more distant among primates.
species_lca = {
    "aotNan1":  {"common_name": "Ma's night monkey", "lca_mya": 33, "rank": 5},
    "calJac3":  {"common_name": "Common marmoset", "lca_mya": 33, "rank": 5},
    "cebCap1":  {"common_name": "Capuchin monkey", "lca_mya": 33, "rank": 5},
    "cerAty1":  {"common_name": "Sooty mangabey", "lca_mya": 33, "rank": 5},
    "chlSab2":  {"common_name": "Green monkey", "lca_mya": 23, "rank": 4},
    "colAng1":  {"common_name": "Angolan colobus", "lca_mya": 23, "rank": 4},
    "gorGor5":  {"common_name": "Gorilla", "lca_mya": 7, "rank": 1},
    "macFas5":  {"common_name": "Crab-eating macaque", "lca_mya": 23, "rank": 4},
    "macNem1":  {"common_name": "Pig-tailed macaque", "lca_mya": 23, "rank": 4},
    "manLeu1":  {"common_name": "Drill", "lca_mya": 23, "rank": 4},
    "nasLar1":  {"common_name": "Proboscis monkey", "lca_mya": 23, "rank": 4},
    "nomLeu3":  {"common_name": "Northern white-cheeked gibbon", "lca_mya": 19, "rank": 3},
    "panPan2":  {"common_name": "Bonobo", "lca_mya": 6, "rank": 1},
    "panTro5":  {"common_name": "Chimpanzee", "lca_mya": 6, "rank": 1},
    "papAnu3":  {"common_name": "Olive baboon", "lca_mya": 23, "rank": 4},
    "ponAbe2":  {"common_name": "Sumatran orangutan", "lca_mya": 13, "rank": 2},
    "rheMac8":  {"common_name": "Rhesus macaque", "lca_mya": 23, "rank": 4},
    "rhiBie1":  {"common_name": "Black snub-nosed monkey", "lca_mya": 23, "rank": 4},
    "rhiRox1":  {"common_name": "Golden snub-nosed monkey", "lca_mya": 23, "rank": 4},
    "saiBol1":  {"common_name": "Bolivian squirrel monkey", "lca_mya": 33, "rank": 5},
}


# ------------------- HELPER FUNCTIONS -------------------
def parse_coord(s):
    """Parse genomic coordinate like 'chr1:1412251.0-1412723.0' → (chrom, start, end)."""
    chrom, rest = s.split(":")
    chrom = chrom.replace("chr", "")
    start_s, end_s = rest.split("-")
    return chrom, int(float(start_s)), int(float(end_s))


def presence_species_for_ltr(ltr_id, ltr_df):
    """Return list of species where given LTR is present above threshold."""
    rows = ltr_df[ltr_df["ltr"] == ltr_id]
    present = rows[rows["coverage"] > COVERAGE_THRESHOLD]["species"].unique().tolist()
    return present


def get_best_species(species_found):
    """Return the most distant species (highest rank) and its MYA."""
    best_species, best_rank = None, -1
    for s in species_found:
        if s in species_lca:
            rank = species_lca[s]["rank"]
            if rank > best_rank:
                best_species, best_rank = s, rank
    best_mya = species_lca[best_species]["lca_mya"] if best_species else None
    return best_species, best_rank, best_mya


# ------------------- MAIN PROCESS -------------------
def main():
    # --- Load data ---
    erv_df = pd.read_csv(
        ERV_BED,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "id", "score", "strand"],
        dtype={"chrom": str},
    )

    paired_df = pd.read_csv(
        PAIRED_LTR,
        sep=r"\s+",
        header=None,
        names=["chr", "ltr1", "ltr2", "len1", "len2", "similarity", "ltr_name"],
    )

    ltr_df = pd.read_csv(LTR_PRESENCE)

    # --- Normalize coordinates ---
    paired_df[["ltr1_chr", "ltr1_start", "ltr1_end"]] = paired_df["ltr1"].apply(
        lambda x: pd.Series(parse_coord(x))
    )
    paired_df[["ltr2_chr", "ltr2_start", "ltr2_end"]] = paired_df["ltr2"].apply(
        lambda x: pd.Series(parse_coord(x))
    )
    erv_df["chrom"] = erv_df["chrom"].astype(str).str.replace("^chr", "", regex=True)

    # --- Iterate through ERVs ---
    results = []

    for _, erv in erv_df.iterrows():
        chrom = str(erv["chrom"])
        if chrom == "Y":  # skip Y chromosome
            continue

        erv_start, erv_end = int(erv["start"]), int(erv["end"])
        window_start, window_end = erv_start - MARGIN, erv_end + MARGIN

        def overlaps(row):
            """Check if an LTR pair overlaps the ERV window."""
            cond1 = (
                row["ltr1_chr"] == chrom
                and not (row["ltr1_end"] < window_start or row["ltr1_start"] > window_end)
            )
            cond2 = (
                row["ltr2_chr"] == chrom
                and not (row["ltr2_end"] < window_start or row["ltr2_start"] > window_end)
            )
            return cond1 or cond2

        matched = paired_df[paired_df.apply(overlaps, axis=1)]
        if matched.empty:
            continue

        species_found = set()
        similarities = []

        for _, p in matched.iterrows():
            ltrs = [p["ltr1"], p["ltr2"]]
            similarities.append(p["similarity"])
            for ltr in ltrs:
                species_found.update(presence_species_for_ltr(ltr, ltr_df))

        # If no non-human species found
        if not species_found or species_found == {"human"}:
            if not similarities:
                continue
            results.append({
                "erv_id": erv["id"],
                "chrom": chrom,
                "erv_start": erv_start,
                "erv_end": erv_end,
                "matched_pairs": len(matched),
                "most_distant_species": "human_only",
                "most_distant_rank": 0,
                "most_distant_mya": 0,
                "similarities": ";".join([f"{s:.3f}" for s in similarities]),
            })
            continue

        best_species, best_rank, best_mya = get_best_species(species_found)
        results.append({
            "erv_id": erv["id"],
            "chrom": chrom,
            "erv_start": erv_start,
            "erv_end": erv_end,
            "matched_pairs": len(matched),
            "most_distant_species": best_species or "unknown",
            "most_distant_rank": best_rank if best_species else 0,
            "most_distant_mya": best_mya if best_species else 0,
            "similarities": ";".join([f"{s:.3f}" for s in similarities]),
        })

    # --- Save results ---
    out_df = pd.DataFrame(results)
    out_df.to_csv(OUTPUT, index=False)
    print(f"Wrote {len(out_df)} records to {OUTPUT}")


if __name__ == "__main__":
    main()

