# ============================================================
# Rarefaction curve comparison across clusters
#
# Analyses:
#   Negative-binomial GAMMs on newly discovered variants
#
# Analysis restricted to ranks 1-500.
# TB_NT is excluded from all analyses.
# ============================================================

# Load libraries
library(tidyverse)
library(mgcv)
library(Cairo)

# Folder containing the Cluster-wise Summed Counts 
data_dir = "../ClusterFiles_06052026"

#Restrict analysis to rank 500
MAX_RANK = 500

#Level of Significance
ALPHA = 0.05

#Exclude Clusters
EXCLUDE_CLUSTERS = c("TB_NT")

# ============================================================
# Read Cluster-wise Input files
# ============================================================
files = list.files(
  data_dir,
  pattern = "_summed_direct_counts(\\..*)?$",
  full.names = TRUE
)

if (length(files) == 0) {
  stop("No files found. Check data_dir and file-name pattern.")
}

# ============================================================
# Function to read one cluster file
# ============================================================
read_cluster_file = function(file) {
  
  cluster_name = basename(file)
  cluster_name = sub("_summed_direct_counts(\\..*)?$", "", cluster_name)
  
  message("Reading: ", basename(file), " | cluster = ", cluster_name)
  
  df = read.table(
    file,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  names(df) = trimws(names(df))
  
  #convert dataframes to long format
  df_long = df %>%
    mutate(cluster = cluster_name, .before = 1) %>%
    pivot_longer(
      cols = -c(cluster, Permutation),
      names_to = "rank",
      values_to = "count"
    ) %>%
    mutate(
      cluster = as.character(cluster),
      iteration = as.integer(Permutation),
      
      # Clean rank names such as rank_1, rank 1, or 1
      rank = as.character(rank),
      rank = sub("^rank[_ ]?", "", rank),
      rank_num = as.integer(rank),
      
      # Cumulative rarefaction count
      count = as.numeric(count)
    )
  
  return(df_long)
}

# ============================================================
# Read all Clusters, combine, exclude cluster(if any) and 
# restrict all analysis till max rank 500 
# ============================================================
combined_df_all = purrr::map_dfr(files, read_cluster_file)

combined_df = combined_df_all %>%
  filter(
    !cluster %in% EXCLUDE_CLUSTERS,
    !is.na(rank_num),
    rank_num <= MAX_RANK
  ) %>%
  mutate(
    cluster = factor(cluster),
    iteration = factor(iteration),
    curve_id = interaction(cluster, iteration, drop = TRUE),
  ) %>%
  droplevels()

if (n_distinct(combined_df$cluster) < 2) {
  stop("Need at least two clusters after exclusion.")
}

message("Clusters retained:")
print(sort(unique(as.character(combined_df$cluster))))

message("Clusters excluded:")
print(EXCLUDE_CLUSTERS)


# ============================================================
# Compute newly discovered variants
# ============================================================
# count        = cumulative rarefaction count (dependent across ranks)
# new_variants = incremental newly discovered variants at each rank (independent)
# ============================================================

combined_df2 = combined_df %>%
  arrange(cluster, iteration, rank_num) %>%
  group_by(curve_id) %>%
  mutate(
    new_variants = count - lag(count, default = 0)
  ) %>%
  ungroup()

if (any(combined_df2$new_variants < 0, na.rm = TRUE)) {
  stop("Negative new_variants detected. Check!!")
}


write.csv(combined_df2,"combined_df2_incremental_counts.csv",row.names = FALSE)



# ============================================================
# Qs: Do the rarefaction curves differ significantly
# ============================================================
# Outcome:
#   new_variants = incremental number of newly discovered variants
#
# Base model: 
#   common smooth rank effect across clusters
#
# Interaction model:
#   cluster-specific smooth rank effects
# ============================================================

m_base_nb = bam(
  new_variants ~ cluster +
    s(rank_num, k = 3) +
    s(curve_id, bs = "re"),
  data = combined_df2,
  family = nb(),
  method = "fREML"
)

m_interaction_nb = bam(
  new_variants ~ cluster +
    s(rank_num, k = 3) +
    s(rank_num, by = cluster, k = 3) +
    s(curve_id, bs = "re"),
  data = combined_df2,
  family = nb(),
  method = "fREML"
)

saveRDS(
  m_base_nb,
  file = "Base_model_Old_k_3.rds"
)

# Save interaction model
saveRDS(
  m_interaction_nb,
  file = "Interaction_model_Old_k_3.rds"
)

# ============================================================
# Model comparison and summaries
# ============================================================
anova_nb = anova(
  m_base_nb,
  m_interaction_nb,
  test = "Chisq"
)

summary_base_nb = summary(m_base_nb)
summary_interaction_nb = summary(m_interaction_nb)

k.check(m_interaction_nb)


# ============================================================
# Residual plots
# ============================================================
diagnostic_df = combined_df2 %>%
  mutate(
    fitted_value = fitted(m_interaction_nb),
    pearson_resid = residuals(m_interaction_nb, type = "pearson"),
    deviance_resid = residuals(m_interaction_nb, type = "deviance")
  )

print("Hello")

#Residual vs sampling rank
p_resid_rank = ggplot(
  diagnostic_df,
  aes(x = rank_num, y = pearson_resid)
) +
  geom_point(alpha = 0.25, size = 0.8) +
  geom_smooth(se = FALSE, linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = "Sampling rank",
    y = "Pearson residual",
    title = "Pearson residuals across sampling rank"
  ) +
  theme_bw()

ggsave(
  "Model_residuals_across_rank.png",
  p_resid_rank,
  width = 7,
  height = 5,
  dpi = 300
)

p_resid_dev_rank = ggplot(
  diagnostic_df,
  aes(x = rank_num, y = deviance_resid)
) +
  geom_point(alpha = 0.25, size = 0.8) +
  geom_smooth(se = FALSE, linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = "Sampling rank",
    y = "Deviance residual",
    title = "Deviance residuals across sampling rank"
  ) +
  theme_bw()

ggsave(
  "Model_residuals_dev_across_rank.png",
  p_resid_dev_rank,
  width = 7,
  height = 5,
  dpi = 300
)


# Histogram of Pearson residuals
p_resid_hist = ggplot(
  diagnostic_df,
  aes(x = pearson_resid)
) +
  geom_histogram(
    bins = 60,
    color = "black",
    fill = "grey70"
  ) +
  labs(
    x = "Pearson residual",
    y = "Frequency",
    title = "Histogram of Pearson residuals"
  ) +
  theme_bw()

ggsave(
  "Model_residuals_histogram.png",
  p_resid_hist,
  width = 7,
  height = 5,
  dpi = 300
)

# Histogram of deviance residuals
p_deviance_hist = ggplot(
  diagnostic_df,
  aes(x = deviance_resid)
) +
  geom_histogram(
    bins = 60,
    color = "black",
    fill = "grey70"
  ) +
  labs(
    x = "Deviance residual",
    y = "Frequency",
    title = "Histogram of deviance residuals"
  ) +
  theme_bw()

ggsave(
  "Model_deviance_residuals_histogram.png",
  p_deviance_hist,
  width = 7,
  height = 5,
  dpi = 300
)

# Residual vs Fitted
p_resid_fitted = ggplot(
  diagnostic_df,
  aes(x = fitted_value, y = deviance_resid)
) +
  geom_point(alpha = 0.25, size = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = "Fitted value",
    y = "Deviance residual",
    title = "Deviance residuals vs fitted values"
  ) +
  theme_bw()

ggsave(
  "Model_deviance_residuals_vs_fitted.png",
  p_resid_fitted,
  width = 6,
  height = 6,
  dpi = 300
)

#=============================================================
# Write Output
#=============================================================
anova_nb
summary_base_nb
summary_interaction_nb

capture.output(
  anova_nb,
  file = "Model_comparison_ranks_1_to_500_TB_NT_excluded.txt"
)

capture.output(
  summary_base_nb,
  file = "Base_Model_summary_ranks_1_to_500_TB_NT_excluded.txt"
)

capture.output(
  summary_interaction_nb,
  file = "Interaction_model_summary_ranks_1_to_500_TB_NT_excluded.txt"
)

message("Residual plots saved.")
message("Analysis complete.")