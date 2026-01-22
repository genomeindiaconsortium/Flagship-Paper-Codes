import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

df = pd.read_excel("PlotInputs.xlsx",sheet_name="Figure1D") 

# Unique group labels
groups = sorted(set(df["Group1"]) | set(df["Group2"]))

# Initialize matrices
jac_mat = pd.DataFrame(np.nan, index=groups, columns=groups)
sample_mat = pd.DataFrame(np.nan, index=groups, columns=groups)

# Fill matrices
for _, row in df.iterrows():
    g1, g2 = row["Group1"], row["Group2"]
    jac = row["Jaccard_values"]
    # samples = row["CombinedSamples"]
    jac_mat.loc[g1, g2] = jac
    jac_mat.loc[g2, g1] = jac  # for symmetry
    # sample_mat.loc[g2, g1] = samples  # upper triangle

# Mask upper triangle for heatmap
mask = np.triu(np.ones_like(jac_mat, dtype=bool))

colors = [
    (0.0, "#ffffcc"),   # light yellow (low values)
    (0.475, "#a1dab4"), # greenish (around 0.475)
    (0.49, "#41b6c4"),  # teal highlight at 0.49
    (0.500, "#2c7fb8"), # blue for 0.500
    (1.0, "#253494")    # dark blue (high values)
]

custom_cmap = LinearSegmentedColormap.from_list("custom", colors)

fig, ax = plt.subplots(figsize=(10, 8))
sns.heatmap(
    jac_mat,
    mask=mask,
    annot=True,
    annot_kws={'size': 15},
    fmt=".3f",
    cmap=custom_cmap,
    square=True,
    linewidths=0.5,
    cbar_kws={"label": "Jaccard Index"},
    ax=ax
)
# Annotate upper triangle manually with sample size
for i in range(len(groups)):
    for j in range(i+1, len(groups)):
        value = sample_mat.iloc[i, j]
        if pd.notna(value):
            ax.text(j + 0.5, i + 0.5, f"{int(value)}", 
                    ha='center', va='center', color='black', fontsize=8)

# Ticks and title
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
plt.tight_layout()
plt.savefig("Figure2D_LowerTriangle.png", dpi=300)
plt.show()
