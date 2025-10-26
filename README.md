
# ERV–LTR Evolutionary Analysis Pipeline

This repository contains a reproducible pipeline for analyzing endogenous retrovirus (ERV) loci, identifying paired LTRs, and quantifying their similarity and evolutionary conservation across primate species.

---

## 📁 Project Structure

```

scripts/
├── download_data.sh       # Download reference genomes and annotations
├── filter_maf.sh          # Extract relevant genomic regions from MAF alignments

src/
├── filter_maf.py          # Parse and filter multiple alignment data
├── find_ltr_pairs.py      # Identify paired LTRs and compute sequence similarity
├── check_erv_presence.py  # Detect presence/absence of LTRs across species
├── summarize_erv_similarity.py  # Summarize ERV–LTR matches and cross-species data
├── plot_erv_similarity_vs_distance.py # Generate plots of LTR similarity vs. evolutionary distance

```

All processed outputs are written to:

```

data/processed/

````
---

## 🚀 Usage

### 1. Download and prepare data

```bash
bash scripts/download_data.sh
bash scripts/filter_maf.sh
```

### 2. Identify paired LTRs

```bash
python src/find_ltr_pairs.py
```

### 3. Check ERV presence across species

```bash
python src/check_erv_presence.py
```

### 4. Summarize ERV–LTR similarity and cross-species presence

```bash
python src/summarize_erv_similarity.py
```

### 5. Plot LTR similarity vs. evolutionary distance

```bash
python src/plot_erv_similarity_vs_distance.py
```

---

## 📊 Outputs

Key output files in `data/processed/`:

- `paired_ltr_results.txt` — nearby LTR pairs and their sequence similarity.
- `ltr_presence.csv` — presence/coverage of each LTR across primate species.
- `erv_ltr_presence_summary.csv` — summary of ERVs with most distant species and LTR similarities.
