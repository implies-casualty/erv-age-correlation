import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Define rank labels
RANK_LABELS = {
    0: 'Human-only',
    1: 'Chimp, Gorilla',
    2: 'Orangutan',
    3: 'Gibbon',
    4: 'Old World Monkeys',
    5: 'New World Monkeys'
}

def parse_args():
    parser = argparse.ArgumentParser('Plot LTR similarity versus evolutionary distance using summarized ERV data.')
    parser.add_argument(
        "--similarity-threshold",
        type=float,
        default=0.8,
        help="Similarity threshold to filter out noise (default: 0.8)",
    ) 
    return parser.parse_args()

args = parse_args()

# Load the data
file_path = 'data/processed/erv_ltr_presence_summary.csv'
data = pd.read_csv(file_path)

# Convert similarities to numeric, coercing errors to NaN
data['similarities'] = pd.to_numeric(data['similarities'], errors='coerce')

# Filter out rows with similarities < threshold
filtered_data = data[data['similarities'] > args.similarity_threshold]

# Group by 'most_distant_rank' and calculate median, average, and count
grouped_stats = filtered_data.groupby('most_distant_rank')['similarities'].agg(['median', 'mean', 'count']).reset_index()

# Keep only ranks with count >= 15
grouped_stats = grouped_stats[grouped_stats['count'] >= 15]

# Function to calculate the confidence interval
def confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.mean(a), np.std(a, ddof=1) / np.sqrt(n)
    h = se * 1.96  # 1.96 is the z-score for 95% confidence interval for large samples
    return m - h, m + h

# Calculate confidence intervals for each rank
confidence_intervals = filtered_data.groupby('most_distant_rank')['similarities'].apply(confidence_interval).apply(pd.Series).reset_index()
confidence_intervals.columns = ['most_distant_rank', 'ci_lower', 'ci_upper']

# Merge the statistics with the confidence intervals
final_stats = pd.merge(grouped_stats, confidence_intervals, on='most_distant_rank')

# --- sort by rank to ensure correct label alignment ---
final_stats = final_stats.sort_values('most_distant_rank').reset_index(drop=True)

# --- Table generation ---

# Define mapping of most_distant_rank -> Last common ancestor
LCA_LABELS = {
    0: "< 6 MYA",           # Human-only
    1: "6–8 MYA",           # Chimp, Gorilla
    2: "12–16 MYA",         # Orangutan
    3: "18–20 MYA",         # Gibbon
    4: "25–30 MYA",         # Old World Monkeys
    5: "35–40 MYA",         # New World Monkeys
}

# Build table DataFrame
table = pd.DataFrame({
    "Most distant group": [RANK_LABELS[rank] for rank in final_stats['most_distant_rank']],
    "Last common ancestor": [LCA_LABELS.get(rank, "?") for rank in final_stats['most_distant_rank']],
    "Average LTR-LTR similarity (95% CI)": [
        f"{row['mean']:.3f} ({row['ci_lower']:.3f}–{row['ci_upper']:.3f})"
        for _, row in final_stats.iterrows()
    ]
})

print("\nGenerated summary table:\n")
print(table.to_markdown(index=False))

# Plotting
plt.figure(figsize=(10, 6))
plt.errorbar(final_stats['most_distant_rank'], final_stats['mean'],
             yerr=[final_stats['mean'] - final_stats['ci_lower'], final_stats['ci_upper'] - final_stats['mean']],
             fmt='-o', ecolor='red', capsize=5) # , linestyle='')

# Adding median and count as text on the plot
for i in range(len(final_stats)):
    plt.text(final_stats['most_distant_rank'][i], final_stats['mean'][i] + 0.01,
             f"Median: {final_stats['median'][i]:.3f}\n#ERVs: {final_stats['count'][i]}",
             ha='center')

plt.xticks(
    ticks=final_stats['most_distant_rank'],
    labels=[RANK_LABELS.get(rank, str(rank)) for rank in final_stats['most_distant_rank']],
    rotation=30,
    ha='right'
)

plt.xlabel('Most Distant Primate')
plt.ylabel('Average Similarity (95% confidence intervals)')
plt.title('Human LTR–LTR Similarity vs. Most Distant Primate Sharing the ERV')
plt.grid(True)
plt.tight_layout()  # Prevent label cutoff
plt.show()
