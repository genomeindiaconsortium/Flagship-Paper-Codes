################################################################################
# Star Allele Metabolizer Phenotype Distribution Across Indian Linguistic Groups
#
# This script generates:
#
#   Figure 4E:
#     Overall metabolizer phenotype distribution across pharmacogenes
#
#   Figure 4F:
#     Population-wise metabolizer phenotype distribution of the
#     most variable genes stratified by linguistic/tribal groups
#
# INPUT DATA:
# Although the original analysis uses Stargazer/Cyrius output files,
# the plots can be reproduced directly from the processed dataframe:
#
#     pop_gene_phenotype (Fig4ef_PGx_star_alleles_table.tsv)
#
# Each row in the dataframe should represent:
#
#     One phenotype category × one gene × one linguistic group
#
# REQUIRED COLUMNS:
#
# 1. Linguistic_group_Tribe
#    - Linguistic/tribal classification of the population.
#    - Allowed values:
#         AA_T
#         DR_NT
#         DR_T
#         IE_NT
#         IE_T
#         TB_NT
#         TB_T
#         CAO
#
# 2. gene
#    - Pharmacogene symbol.
#    - Should preferably be uppercase.
#    - Example:
#         CYP2D6
#         CYP2C19
#         SLCO1B1
#
# 3. phenotype
#    - Functional/metabolizer phenotype category.
#    - Allowed values after harmonization:
#         Normal
#         Intermediate
#         Poor
#         Rapid
#         Unknown
#
# 4. count
#    - Number of individuals belonging to the corresponding phenotype
#      category within that linguistic group.
#
# 5. freq
#    - Percentage frequency of the phenotype within that linguistic group.
#    - Calculated as:
#
#         freq = (count / total individuals for that gene and group) × 100
#
# DATA PROCESSING OVERVIEW:
#
# 1. Stargazer genotype outputs were parsed for all pharmacogenes.
#
# 2. CYP2D6 phenotypes generated using Cyrius were merged separately
#    due to improved structural variant calling accuracy.
#
# 3. Individual-level phenotypes were merged with population metadata
#    containing linguistic/tribal classifications.
#
# 4. Raw phenotype labels from Stargazer/Cyrius were harmonized into
#    unified phenotype classes:
#
#         Normal
#         Intermediate
#         Poor
#         Rapid
#         Unknown
#
# 5. Phenotype frequencies were calculated for each:
#
#         Linguistic group × Gene combination
#
# FIGURE OUTPUTS:
#
# Figure 4E (object: p)
# -----------------------------------------
# Horizontal stacked barplot showing overall phenotype distribution
# across pharmacogenes.
#
# Figure 4F (object: q)
# -----------------------------------------
# Faceted stacked barplot showing phenotype distributions across
# linguistic/tribal groups for selected most variable genes.
#
# SELECTED GENES IN FIGURE 4F:
#
#     CYP2B6
#     CYP2C19
#     CYP2D6
#     CYP3A5
#     CYP4F2
#     SLCO1B1
#     UGT1A1
#     VKORC1
#
# OUTPUT FILES:
#
#     gene_matrix_main.png        -> Figure 4E
#     pop_gene_matrix_main.png    -> Figure 4F
#
################################################################################

library(ggplot2)
library(tidyverse)
library(readr)
library(readxl)

cyp_genes <- c("cyp2b6", "cyp2c9", "cyp3a5", "cyp2c19", "cyp4f2", "cyp2d6")
non_cyp_genes <- c("cftr", "slco1b1", "dpyd", "tpmt", "ifnl3", "ugt1a1", "nudt15", "vkorc1")

pop_info <- read_excel("final_pop_info.xlsx") %>%
  select(FID, Linguistic_group_Tribe) %>% 
  distinct()

cyp2d6 <- read_excel("cyp2d6_merged_genotypes.xlsx") %>%
  rename(name = FID, phenotype = Metabolizer_Status) %>% 
  select(name, Linguistic_group_Tribe, phenotype) %>% 
  mutate(gene = "cyp2d6")

cyrius_samples <- cyp2d6$name

merged_rows_df <- data.frame()

for (gene in c(cyp_genes, non_cyp_genes)) {
  star_allele <- file.path("stargazer/", paste0("output-", gene, ".sg-genotype.txt"))
  
  if (!file.exists(star_allele)) {
    message("Gene not genotyped: ", gene)
    next
  }
  
  df <- read_tsv(star_allele, col_types = cols(.default = "c")) %>%
    select(name, phenotype) %>% 
    mutate(gene = gene)
  
  if(gene == "cyp2d6") {
  df <- df %>%
    filter(!name %in% cyrius_samples)
  }
  
  merged_df <- merge(df, pop_info, by.x = "name", by.y = "FID") %>%
    filter(Linguistic_group_Tribe != "")
  
  merged_rows_df <- bind_rows(merged_rows_df, merged_df)
}

final_df <- bind_rows(merged_rows_df, cyp2d6)

pop_gene_phenotype <- final_df %>%
  group_by(Linguistic_group_Tribe, gene, phenotype) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Linguistic_group_Tribe, gene) %>%
  mutate(freq = (count / sum(count)) * 100)

pop_gene_phenotype <- pop_gene_phenotype %>%
  mutate(phenotype = recode(phenotype,
                            "normal_metabolizer" = "Normal",
                            "Normal Metabolizer" = "Normal",
                            "normal_function" = "Normal",
                            "favorable_response" = "Normal",
                            "intermediate_metabolizer" = "Intermediate",
                            "Intermediate Metabolizer" = "Intermediate",
                            "poor_metabolizer" = "Poor",
                            "Poor Metabolizer" = "Poor",
                            "poor_function" = "Poor",
                            "decreased_function" = "Poor",
                            "rapid_metabolizer" = "Rapid",
                            "increased_function" = "Rapid",
                            "ultrarapid_metabolizer" = "Rapid",
                            "Ultrarapid Metabolizer" = "Rapid",
                            "unknown_metabolizer" = "Unknown",
                            "Unknown Metabolizer" = "Unknown",
                            "unknown_function" = "Unknown",
                            "Indeterminate" = "Unknown"),
         Linguistic_group_Tribe = recode(Linguistic_group_Tribe, 
                            "AA_Tribal" = "AA_T",
                            "DR_Non_Tribal" = "DR_NT",
                            "DR_Tribal" = "DR_T",
                            "IE_Non_Tribal" = "IE_NT",
                            "IE_Tribal" = "IE_T",
                            "TB_Non_Tribal" = "TB_NT",
                            "TB_Tribal" = "TB_T",
                            "CAO" = "CAO"),
          Linguistic_group_Tribe = factor(Linguistic_group_Tribe, levels = rev(c("AA_T", "DR_NT", "DR_T", "IE_NT",
                                                                            "IE_T", "TB_NT", "TB_T", "CAO"))),
         gene = toupper(as.character(gene)))

gene_order <- c(
  "CYP2B6",
  "CYP2C9",
  "CYP3A5",
  "CYP2C19",
  "CYP4F2",
  "CYP2D6",
  "CFTR",
  "SLCO1B1",
  "DPYD",
  "TPMT",
  "IFNL3",
  "UGT1A1",
  "NUDT15",
  "VKORC1"
)

pop_gene_phenotype$gene <- factor(
  pop_gene_phenotype$gene,
  levels = gene_order
)

phenotype_colors <- c(
  "Normal" = "#1f77b4",
  "Intermediate" = "#ff7f0e",
  "Poor" = "#d62728",
  "Rapid" = "#2ca02c",
  "Unknown" = "grey70"
)

p <- ggplot(pop_gene_phenotype, aes(x = freq, y = forcats::fct_rev(gene), fill = phenotype)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = phenotype_colors) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Proportion", y = NULL, fill = "Metabolizer") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(face = "bold.italic", size = 6),
    axis.title.x = element_text(size = 8),
    legend.position = "right",
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.3, "cm"),
    legend.title = element_text(size = 6),
    panel.spacing.x = unit(0.08, "cm"),
    panel.spacing.y = unit(-0.02, "cm")
  )

ggsave(
  filename = "gene_matrix_main.png",
  plot = q,
  device = "png",
  width = 3,
  height = 2.8,
  units = "in",
  dpi = 600,
  bg = "white"
)

pop_gene_phenotype <- pop_gene_phenotype %>% filter(gene %in% c("CYP2B6", "CYP2C19", "CYP2D6", "CYP3A5", "CYP4F2", "SLCO1B1", "UGT1A1", "VKORC1"))

q <- ggplot(pop_gene_phenotype, aes(x = freq, y = Linguistic_group_Tribe, fill = phenotype)) +
  geom_bar(stat = "identity", width = 0.8) +
  facet_grid(. ~ gene, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = phenotype_colors) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Proportion", y = NULL, fill = "Metabolizer") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(face = "bold", size = 6),
    axis.title.x = element_text(size = 8),
    strip.background = element_rect(fill = "grey85", color = "grey50"),
    strip.text.x = element_text(face = "bold.italic", size = 6),
    legend.position = "bottom",
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.3, "cm"),
    legend.title = element_text(size = 6),
    panel.spacing.x = unit(0.08, "cm"),
    panel.spacing.y = unit(-0.02, "cm")
  )

ggsave(
  filename = "pop_gene_matrix_main.png",
  plot = r,
  device = "png",
  width = 9,
  height = 2.8,
  units = "in",
  dpi = 600,
  bg = "white"
)