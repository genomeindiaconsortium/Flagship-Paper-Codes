import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.patches import Rectangle
import pandas as pd


# ------------------------------------------------------------------
# Expected structure of merged_df before plotting
#
# merged_df should contain the following columns:
#
#   Population : Population / group name
#   Dataset    : Variant category (overall, dmm, lof)
#   Common     : Number of common variants
#   Rare       : Number of rare variants
#
# Example:
#
#   Population    Dataset    Common      Rare
#   ------------------------------------------------
#   DR_Tribe      overall    9218619     10170956
#   DR_NonTribe   overall    8941521     25034224
#   IE_NonTribe   overall    8931914     45366127
#   ...
#
# - Each population should have one row per dataset category.
# - Dataset categories used in this plot are:
#       overall
#       dmm
#       lof
# - Common and Rare columns must contain numeric variant counts.
# - Population names will later be shortened using rename_map.
#
# Final merged_df format should look like:
#
#       Population      Dataset     Common      Rare
#   0   DR_Tribe        overall     9218619     10170956
#   1   DR_NonTribe     overall     8941521     25034224
#   2   IE_NonTribe     overall     8931914     45366127
#   3   IE_Tribe        overall     8949734     17607735
#   4   AA_Tribe        overall     8849430     13370002
#   5   TB_NonTribe     overall     8403329     5749457
#   6   TB_Tribe        overall     8490896     12963697
#   7   DR_Tribe        dmm         431         3123
#   8   DR_NonTribe     dmm         322         8998
#   ...
# ------------------------------------------------------------------


prop_df = merged_df.copy()


rename_map = {
    "DR_Tribe":       "DR_T",
    "DR_NonTribe":    "DR_NT",
    "IE_Tribe":       "IE_T",
    "IE_NonTribe":    "IE_NT",
    "AA_Tribe":       "AA_T",
    "TB_Tribe":       "TB_T",
    "TB_NonTribe":    "TB_NT",
}
prop_df["Population"] = prop_df["Population"].replace(rename_map)


GI_add = [
    ["GI", "overall", 8877722, 61869089],
    ["GI", "dmm",      272,    24852],
    ["GI", "lof",      399,    15450]
]
prop_df = pd.concat([prop_df, pd.DataFrame(GI_add, columns=["Population","Dataset","Common","Rare"])], ignore_index=True)


prop_df["Total"] = prop_df["Common"] + prop_df["Rare"]
prop_df["Common_prop"] = prop_df["Common"] / prop_df["Total"]
prop_df["Rare_prop"]   = prop_df["Rare"] / prop_df["Total"]

datasets = ["overall", "dmm", "lof"]
dataset_labels = {"overall": "Overall", "dmm": "DMM", "lof": "LOF"}


sorted_populations = ["GI", "AA_T", "DR_NT", "DR_T", "IE_NT", "IE_T", "TB_NT", "TB_T"]


color_map = {
    "overall": ("#BDBDBD", "#6EB1E6"),
    "dmm":     ("#969696", "#2589DA"),
    "lof":     ("#737373", "#0E4675")
}


sns.set_theme(style="white", context="talk")
plt.rcParams.update({
    'font.size': 18,
    'axes.labelsize': 18,
    'axes.titlesize': 20,
    'xtick.labelsize': 16,
    'ytick.labelsize': 18,
    'legend.fontsize': 16
})

fig, ax = plt.subplots(figsize=(9, 8))

y = np.arange(len(sorted_populations))
height = 0.25


for i, dataset in enumerate(datasets):
    subset = prop_df[prop_df["Dataset"] == dataset].set_index("Population").reindex(sorted_populations)

    rare_vals   = subset["Rare_prop"].fillna(0)
    common_vals = subset["Common_prop"].fillna(0)
    total_vals  = subset["Total"].fillna(0)

    bar_y = y + (i - 1) * height
    common_color, rare_color = color_map[dataset]

    bars_rare = ax.barh(
        bar_y, rare_vals, height,
        color=rare_color,
        edgecolor='white', linewidth=1.0
    )

    bars_common = ax.barh(
        bar_y, common_vals, height,
        left=rare_vals,
        color=common_color,
        edgecolor='white', linewidth=1.0
    )

    
    for rect, val in zip(bars_rare, rare_vals):
        if val > 0.06:
            ax.text(
                rect.get_width()/2,
                rect.get_y() + rect.get_height()/2,
                f"{val:.1%}",
                ha='center', va='center',
                color='white', fontsize=14, fontweight='bold'
            )

    
    for j, (rv, cv, y_pos) in enumerate(zip(rare_vals, common_vals, bar_y)):
        ax.text(
            1.01,                           # slightly beyond bar end (proportion scale)
            y_pos,
            f"{int(total_vals.iloc[j]):,}",
            va='center',
            ha='left',
            fontsize=13,
            fontweight='bold',
            color='black'
        )


ax.set_yticks(y)
ax.set_yticklabels(sorted_populations, fontsize=18, fontweight='bold')
ax.set_xlabel("Proportion", fontweight='bold')
ax.set_xlim(0, 1.12)   # space for totals
ax.invert_yaxis()      # <- this places the first entry (GI) at the top

sns.despine(left=True)
ax.xaxis.grid(True, linestyle='--', alpha=0.5)
ax.set_axisbelow(True)


legend_elements = []
for ds in datasets:
    cc, rc = color_map[ds]
    legend_elements.append(Rectangle((0,0),1,1, color=rc, label=f"{dataset_labels[ds]} (Rare)"))
    legend_elements.append(Rectangle((0,0),1,1, color=cc, label=f"{dataset_labels[ds]} (Common)"))

ax.legend(
    handles=legend_elements,
    loc='lower center',
    bbox_to_anchor=(0.5, -0.35),
    ncol=3,
    frameon=True,
    facecolor='white',
    edgecolor='lightgray',
    fontsize=16
)

plt.tight_layout()
plt.savefig("Horizontal_bar_4A_final_GI.png", dpi=500, bbox_inches='tight')
plt.show()