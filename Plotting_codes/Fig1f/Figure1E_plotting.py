import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

input_dir = "ClusterFiles" 

'''
Input directory contains cluster-wise counts - file name: <cluster>_summed_direct_counts.tsv 

An example file for AA_T_summed_direct_counts.tsv look like this:

#===============
Iteration	rank_1	rank_2	rank_3	...   rank_500

1	58861	80957	94786	104892  ...    1462911
2	58107	77936	92327	101746	...    1461390
3	58167	79125	91964	102608	...    1463970
...
...

100	59476	82001	95636	107859  ...    1468733
#================

'''

output_png = "Figure1F_withoutLegend.png"

LOWER_Q = 2.5
UPPER_Q = 97.5
MAX_POINTS = 500

EXCLUDE_CLUSTERS = {"TB_NT"}

cluster_colors = {
    "IE_T":  "#cd3333",
    "DR_T":  "#000080",
    "AA_T":  "#4682b4",
    "TB_T":  "#7ccd7c",
    "IE_NT": "#FAC4C4",
    "DR_NT": "#8A8AF4",
    "TB_NT": "#C7F2C7"
}


TITLE_FSIZE = 24
LABEL_FSIZE = 22
TICK_FSIZE = 20
LEGEND_FSIZE = 20


# For each cluster, this script reads the cumulative rarefaction counts
# across 100 rarefaction iterations. 
# At each sampling rank, it computes the mean cumulative count and 
# the 2.5th-97.5th percentile interval across replicates representing variability due to rarefaction resampling and 
# then plots these cluster-specific rarefaction curves.

def summarize_summed_counts(file_path, lower_q=2.5, upper_q=97.5, max_points=None):
    df = pd.read_csv(file_path, sep="\t")

    if "Permutation" not in df.columns:
        raise ValueError(f"'Permutation' column missing in {file_path}")

    sample_cols = [c for c in df.columns if c != "Permutation"]

    if max_points is not None:
        sample_cols = sample_cols[:max_points]

    rows = []

    for idx, sample_name in enumerate(sample_cols, start=1):
        vals = pd.to_numeric(df[sample_name], errors="coerce").dropna().values

        if len(vals) == 0:
            continue

        rows.append({
            "rank": idx,
            "sample_name": sample_name,
            "mean": np.mean(vals),
            "ci_low": np.percentile(vals, lower_q),
            "ci_high": np.percentile(vals, upper_q)
        })

    return pd.DataFrame(rows)


# The code will look for files ending with '_summed_direct_counts.tsv' in input directory
summed_files = sorted([
    os.path.join(input_dir, f)
    for f in os.listdir(input_dir)
    if f.endswith("_summed_direct_counts.tsv")
])

if not summed_files:
    raise ValueError(f"No *_summed_direct_counts.tsv files found in {input_dir}")

# Plot figure 
plt.figure(figsize=(14, 8))

for file_path in summed_files:
    fname = os.path.basename(file_path)
    cluster_name = re.sub(r"_summed_direct_counts\.tsv$", "", fname)

    if cluster_name in EXCLUDE_CLUSTERS:
        print(f"Skipping {cluster_name}")
        continue

    summary_df = summarize_summed_counts(
        file_path,
        lower_q=LOWER_Q,
        upper_q=UPPER_Q,
        max_points=MAX_POINTS
    )

    color = cluster_colors.get(cluster_name, "#808080")

    plt.fill_between(
        summary_df["rank"],
        summary_df["ci_low"],
        summary_df["ci_high"],
        color=color,
        alpha=0.30
    )

    plt.plot(
        summary_df["rank"],
        summary_df["mean"],
        color=color,
        linewidth=3,
        label=cluster_name
    )

plt.xlabel("Number of samples", fontsize=LABEL_FSIZE)
plt.ylabel("Cumulative novel variants discovered", fontsize=LABEL_FSIZE)
#plt.title("Combined cumulative discovery curves", fontsize=TITLE_FSIZE)

plt.xticks(fontsize=TICK_FSIZE)
plt.yticks(
    [0, 500000, 1000000, 1500000],
    ["0.0M", "0.5M", "1.0M", "1.5M"],
    fontsize=TICK_FSIZE
)


# plt.legend(fontsize=LEGEND_FSIZE)
plt.tight_layout()
plt.savefig(output_png, dpi=300, bbox_inches="tight")
plt.show()

print(f"Saved combined plot to: {output_png}")