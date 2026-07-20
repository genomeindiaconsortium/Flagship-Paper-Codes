Rarefaction_analysis_code.py:
This script takes an input directory containing chromosome-wise summary files (*summaryData_novel.csv). These files provide summary information for each novel variant, including the clusters in which the variant is observed and the corresponding sample names.
The script uses the sample-to-cluster mapping file to identify which sample belongs to which cluster (sample_cluster.tsv) and then performs cluster-wise rarefaction analysis.

The script generates two main outputs for each cluster and chromosome:
1. *_sample_names.tsv
   - Contains the sample order used in each permutation.
   - This file is used to verify that the same sample order is maintained across chromosomes.
2. *_direct_counts.tsv
   - Contains the cumulative novel variant counts for each permutation and rank.
   - Rows represent permutations, and columns represent sample ranks.
This script requires MPI to parallelize processing across multiple chromosomes. The required Python modules are pandas, os, random, collections, and mpi4py.
After running this script, chromosome-wise outputs are summed across chr1–chr22 to generate genome-wide cumulative novel variant counts for each cluster.

Plotting code: The plotting script is used to generate Figure 1E. It takes the merged chromosome-summed files from the "ClusterFiles": directory as input and generates cumulative rarefaction plots.