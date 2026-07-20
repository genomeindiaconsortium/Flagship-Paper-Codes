import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np



# ------------------------------------------------------------------
# Expected structure of df_final before generating the heatmap
#
# df_final should contain:
#
#   Chr                  : Chromosome
#   POS                  : Genomic position
#   gnomADg_AF           : Global gnomAD allele frequency
#   gnomADg_SAS_AF       : gnomAD South Asian allele frequency
#
# Population-specific MAF columns:
#
#   MAF_DR_Tribe_maf.frq
#   MAF_DR_NonTribe_maf.frq
#   MAF_IE_NonTribe_maf.frq
#   MAF_TB_Tribe_maf.frq
#   MAF_IE_Tribe_maf.frq
#   MAF_AA_Tribe_maf.frq
#   MAF_TB_NonTribe_maf.frq
#
# Each row represents one variant.
#
# Example:
#
#   Chr   POS     gnomADg_SAS_AF   MAF_DR_Tribe_maf.frq   ...
#   ----------------------------------------------------------------
#   chr1  963994  0.000207         0.000000
#   chr1  964371  0.000208         0.000000
#   chr1  964409  0.000000         0.000000
#   chr1  964460  0.000414         0.001866
#   ...
#
# Important preprocessing assumptions:
#
# 1. Missing MAF values are replaced with 0:
#
#       df_final[pop_cols] = df_final[pop_cols].fillna(0)
#       df_final['gnomADg_SAS_AF'] = (
#           df_final['gnomADg_SAS_AF'].fillna(0)
#       )
#
# 2. The heatmap visualizes the mean absolute MAF difference
#    between each Indian population and gnomAD SAS.
#
# 3. For every variant:
#
#       abs(population_MAF - gnomADg_SAS_AF)
#
#    is calculated.
#
# 4. Variants are grouped into quantile-based bins according to
#    the maximum observed MAF difference across populations:
#
#       df_final['max_diff']
#
# 5. The final heatmap matrix:
#
#       rows    -> populations
#       columns -> variant bins
#       values  -> mean absolute MAF difference
#
# 6. Population ordering used in the figure:
#
#       DR_Tribe
#       TB_Tribe
#       IE_Tribe
#       AA_Tribe
#       DR_NonTribe
#       IE_NonTribe
#       TB_NonTribe
#
# Final df_final format should look like:
#
#       Chr   POS     gnomADg_SAS_AF   MAF_DR_Tribe_maf.frq
#   0   chr1  963994  0.000207         0.000000
#   1   chr1  964371  0.000208         0.000000
#   2   chr1  964409  0.000000         0.000000
#   3   chr1  964460  0.000414         0.001866
#   4   chr1  964527  0.000000         0.000000
#   ...
# ------------------------------------------------------------------

# Population MAF columns
pop_cols = [
    'MAF_DR_Tribe_maf.frq',
    'MAF_TB_Tribe_maf.frq',
    'MAF_IE_Tribe_maf.frq',
    'MAF_DR_NonTribe_maf.frq',
    'MAF_IE_NonTribe_maf.frq',
    'MAF_TB_NonTribe_maf.frq',
    'MAF_AA_Tribe_maf.frq'  
]

# Fill missing values
df_final[pop_cols] = df_final[pop_cols].fillna(0)
df_final['gnomADg_SAS_AF'] = df_final['gnomADg_SAS_AF'].fillna(0)

# Number of bins
n_bins = 100

# Compute max difference per variant
df_final['max_diff'] = df_final[pop_cols].subtract(df_final['gnomADg_SAS_AF'], axis=0).abs().max(axis=1)

# Create bins
df_final['bin'] = pd.qcut(df_final['max_diff'], n_bins, labels=False, duplicates='drop')

# Aggregate mean MAF difference per bin and population
agg_matrix = df_final.groupby('bin')[pop_cols].apply(lambda x: (x.subtract(df_final.loc[x.index, 'gnomADg_SAS_AF'], axis=0)).abs().mean())

# Reorder rows: populations by Tribe -> Non-Tribe
ordered_pops = [
    'MAF_DR_Tribe_maf.frq',
    'MAF_TB_Tribe_maf.frq',
    'MAF_IE_Tribe_maf.frq',
    'MAF_AA_Tribe_maf.frq',
    'MAF_DR_NonTribe_maf.frq',
    'MAF_IE_NonTribe_maf.frq',
    'MAF_TB_NonTribe_maf.frq'
]

agg_matrix = agg_matrix[ordered_pops].T  # transpose so populations are rows

# Plot
plt.figure(figsize=(8, 6))
sns.heatmap(
    agg_matrix,
    cmap="viridis",
    linewidths=0.5,
    linecolor='gray',
    cbar_kws={'label': 'Mean MAF difference vs gnomAD'},
    vmin=0,
    vmax=np.percentile(agg_matrix.values.flatten(), 95),
    xticklabels=False,
    yticklabels=[pop.replace('MAF_', '').replace('_maf.frq','') for pop in agg_matrix.index]
)

plt.title("Variant-level MAF Differences with gnomAD", fontsize=14, pad=12)
plt.xlabel("Variant bins", fontsize=12)
plt.ylabel("Populations", fontsize=12)
plt.tight_layout()
plt.savefig("Heatmap_DMM_MAF_Diff_gnomAD.png", dpi=600, bbox_inches='tight')
plt.show()