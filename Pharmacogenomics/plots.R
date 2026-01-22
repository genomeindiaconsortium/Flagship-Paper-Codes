library(ggplot2)
library(tidyverse)
library(readr)
library(readxl)
library(ComplexHeatmap)
library(circlize)

# Figure 4d
# Similar for plots S13.1, S13.2, S13.3, S13.4
imp_var <- read_excel("final_pgx_data.xlsx")

imp_var <- imp_var |>
  mutate(populations = as.character(populations)) |>
  mutate(pop_num = suppressWarnings(as.numeric(populations))) |>
  arrange(is.na(pop_num), pop_num) |>
  mutate(populations = factor(populations, levels = c(as.character(1:82), "SAS", "AFR", "AMR", "EAS", "EUR")),
         type = factor(type, levels = c("Toxicity", "Dosage", "Efficacy", "Other")),
         variant_gene = paste0(variant, " (", gene, ")"),
         chemicals = str_to_title(chemicals),
         drug_category = paste0(chemicals, " (", category, ")"),
         Linguistic_group_Tribe = factor(Linguistic_group_Tribe,
                                         levels = c("AA_T", "DR_T", "DR_NT",
                                                    "IE_T", "IE_NT", "TB_T",
                                                    "NT", "Global"))) |>
  arrange(type, Testing, variant_gene) |> 
  group_by(variant_gene) |>
  mutate(norm_frequency = as.numeric(scale(frequency))) |>
  ungroup()

mat <- imp_var |>
  select(populations, variant_gene, norm_frequency) |>
  pivot_wider(names_from = populations, values_from = norm_frequency) |>
  column_to_rownames("variant_gene") |>
  as.matrix()

testing_colors <- c(
  "Testing Required" = "#B30568",
  "Testing Recommended" = "#0568B3",
  "Actionable PGx" = "#69B305"
)
linguistic_colors <- c(
  "AA_T" = "#4682b4", "DR_T" = "#000080", "DR_NT" = "#8A8AF4",
  "IE_T" = "#cd3333", "IE_NT" = "#FAC4C4", "TB_T" = "#7ccd7c",
  "NT" = "#C7F2C7", "Global" = "#000000"
)
type_colors <- c("Toxicity"="#FF9999", "Dosage"="#99CCFF", "Efficacy"="#99FF99", "Other"="#CCCCCC")
bio_geo_colors <- c("Northern Riverine Plains" = "#289D3D", "Western Plains" = "#735323", 
                    "Western Himalayas" = "#C0F1FF", "Western Coastal Plains" = "#9AF708", 
                    "Eastern Riverine Plains" = "#00FFC1", "Eastern Coastal Plains" = "#43582F", 
                    "North Deccan" = "#F8C1A6", "North Central Highlands" = "#EDF181", 
                    "South Central Highlands" = "#A5506D", "South Deccan" = "#BCCB34", 
                    "Eastern Ghats" = "#D1AAC2", "Western Ghats" = "#DB7003", 
                    "Nilgiri Hills" = "#7B9EFF", "Eastern Plateau" = "#FBA600", 
                    "Brahmaputra Valley" = "#328C97", "Central Himalayas" = "#FCFF00", 
                    "North Eastern Range" = "#BBA980", "Global" = "#000000")

min_val <- min(imp_var$norm_frequency, na.rm = TRUE)
max_val <- max(imp_var$norm_frequency, na.rm = TRUE)

col_fun <- colorRamp2(
  c(min_val, -2.5, 0, 2.5, 5, max_val),
  c("#2E627AFF", "#1283c0ff", "#F5F5F5FF", "#D6604DFF", "#B2182BFF", "#67001FFF")
)

bio_geo_data <- imp_var |> distinct(populations, Biogeography)
bio_geo_data <- bio_geo_data[match(colnames(mat), bio_geo_data$populations), ]

row_meta <- imp_var |>
  distinct(variant_gene, drug_category, type, Testing)
row_meta <- row_meta[match(rownames(mat), row_meta$variant_gene), ]

col_meta <- imp_var |>
  distinct(populations, Linguistic_group_Tribe, Biogeography)
col_meta <- col_meta[match(colnames(mat), col_meta$populations), ]

col_ha <- HeatmapAnnotation(
  Linguistic_Group = col_meta$Linguistic_group_Tribe,
  col = list(Linguistic_Group = linguistic_colors),
  annotation_name_gp = gpar(col = NA),
  show_legend = FALSE
)

right_ha <- rowAnnotation(
  Drug = anno_text(row_meta$drug_category,
                   gp = gpar(fontsize = 28, fontface = "italic"),
                   just = "left",
                   location = 0.01)
)

left_ha <- rowAnnotation(
  Type = row_meta$type,
  Variant = anno_text(rownames(mat), 
                      gp = gpar(fontsize = 32, fontface = "bold",
                                col = testing_colors[row_meta$Testing]),
                      just = "right",
                      location = 0.99),
  col = list(Type = type_colors),
  gap = unit(4, "mm"), 
  show_annotation_name = FALSE,
  show_legend = FALSE
)

ht <- Heatmap(
  mat,
  name = "Z-score Frequency",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = TRUE,
  row_split = row_meta$type,
  row_title_gp = gpar(fontsize = 22, fontface = "bold"),
  column_split = col_meta$Linguistic_group_Tribe,
  column_title_gp = gpar(fontsize = 28, fontface = "bold"),
  column_names_gp = gpar(
    fontsize = 30, fontface = "bold",
    col = bio_geo_colors[col_meta$Biogeography]
  ),
  row_gap = unit(2, "mm"),
  column_gap = unit(2, "mm"),
  top_annotation = col_ha,
  right_annotation = right_ha,
  left_annotation = left_ha,
  border = TRUE,
  rect_gp = gpar(col = "white", lwd = 1),
  show_heatmap_legend = FALSE
)

lg_z <- Legend(
  col_fun = col_fun,
  title = "Z-score Frequency",
  at = seq(floor(min_val), ceiling(max_val), by = 2),
  labels = seq(floor(min_val), ceiling(max_val), by = 2),
  direction = "horizontal",
  title_gp = gpar(fontsize = 28, fontface = "bold"),
  labels_gp = gpar(fontsize = 22),
  legend_width = unit(16, "cm"),
  grid_height = unit(8, "mm"),
  title_gap = unit(5, "mm")
)

lg_testing <- Legend(
  labels = names(testing_colors),
  title = "Testing",
  legend_gp = gpar(fill = testing_colors),
  direction = "horizontal",
  title_gp = gpar(fontface = "bold", fontsize = 28),
  labels_gp = gpar(fontsize = 22),
  title_gap = unit(5, "mm"),
  row_gap = unit(2, "mm")
)

combined_legends <- packLegend(
  lg_z, lg_testing,
  direction = "horizontal",
  gap = unit(6, "cm")
)

png("final_heatmap.png",
    width = 60, height = 15, units = "in", res = 600)

draw(
  ht,
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom",
  merge_legend = FALSE,
  heatmap_legend_list = list(combined_legends)
)

dev.off()

# Figure 4e, 4f
cyp_genes <- c("cyp2b6", "cyp2c9", "cyp3a5", "cyp2c19", "cyp4f2")
non_cyp_genes <- c("cftr", "slco1b1", "dpyd", "tpmt", "ifnl3", "ugt1a1", "nudt15", "vkorc1")

pop_info <- read_excel("pop_info.xlsx") |>
  select(FID, Linguistic_group_Tribe) |> 
  distinct()

cyp2d6 <- read_excel("merged_genotypes.xlsx") |>
  rename(name = FID, phenotype = Metabolizer_Status) |> 
  select(name, Linguistic_group_Tribe, phenotype) |> 
  mutate(gene = "cyp2d6")

merged_rows_df <- data.frame()

for (gene in c(cyp_genes, non_cyp_genes)) {
  star_allele <- file.path("stargazer\\",
                           paste0("output-", gene, ".sg-genotype.txt"))
  
  if (!file.exists(star_allele)) {
    message("Gene not genotyped: ", gene)
    next
  }
  
  df <- read_tsv(star_allele, col_types = cols(.default = "c")) |>
    select(name, phenotype) |> 
    mutate(gene = gene)
  
  merged_df <- merge(df, pop_info, by.x = "name", by.y = "FID") |>
    filter(Linguistic_group_Tribe != "")
  
  merged_rows_df <- bind_rows(merged_rows_df, merged_df)
}

final_df <- bind_rows(merged_rows_df, cyp2d6)

pop_gene_phenotype <- final_df |>
  group_by(Linguistic_group_Tribe, gene, phenotype) |>
  summarise(count = n(), .groups = "drop") |>
  group_by(Linguistic_group_Tribe, gene) |>
  mutate(freq = (count / sum(count)) * 100)

pop_gene_phenotype <- pop_gene_phenotype |>
  mutate(phenotype = recode(phenotype,
                            "normal_metabolizer" = "Normal Metabolizer",
                            "normal_function" = "Normal Metabolizer",
                            "favorable_response" = "Normal Metabolizer",
                            "intermediate_metabolizer" = "Intermediate Metabolizer",
                            "poor_metabolizer" = "Poor Metabolizer",
                            "poor_function" = "Poor Metabolizer",
                            "decreased_function" = "Poor Metabolizer",
                            "rapid_metabolizer" = "Rapid Metabolizer",
                            "increased_function" = "Rapid Metabolizer",
                            "ultrarapid_metabolizer" = "Ultrarapid Metabolizer",
                            "unknown_metabolizer" = "Unknown Metabolizer",
                            "unknown_function" = "Unknown Metabolizer",
                            "Indeterminate" = "Unknown Metabolizer"),
         Linguistic_group_Tribe = factor(Linguistic_group_Tribe, levels = c("AA_T", "DR_T", "DR_NT",
                                                                            "IE_T", "IE_NT", "TB_T",
                                                                            "NT")),
         gene = toupper(as.character(gene)))

phenotype_colors <- c(
  "Normal Metabolizer" = "#1f77b4",
  "Intermediate Metabolizer" = "#ff7f0e",
  "Poor Metabolizer" = "#d62728",
  "Rapid Metabolizer" = "#2ca02c",
  "Ultrarapid Metabolizer" = "#9467bd",
  "Unknown Metabolizer" = "grey70"
)

p <- ggplot(pop_gene_phenotype, aes(x = freq, y = gene, fill = phenotype)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(gene ~ Linguistic_group_Tribe, scales = "free_y", space = "free_y") +
  scale_y_discrete(expand = expansion(add = c(0, 0))) +
  scale_fill_manual(values = phenotype_colors) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y.left = element_text(vjust = 0.5, face = "bold", size = 10),
    strip.text.x.top = element_text(face = "bold", size = 10),
    strip.text.y.right = element_blank(),
    legend.position = "right",
    legend.title = element_blank(),
    panel.spacing.x = unit(-0.1, "cm")
  )

ggsave("pop_gene_matrix.jpeg",
       p, device = "jpeg", height = 5, width = 12, dpi = 600)

# Figure S13.4
actionable <- read_excel("actionable_per_sample.xlsx", sheet = "Total_Actionable")
actionable <- actionable |> filter(Population != "Siddi")

# Part A

(each_sample <- ggplot(actionable, aes(x = as.factor(total_actionable))) +
  geom_bar(fill = "#8da0cb", color = "black") +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.3, size = 5) +
  labs(x = "Number of actionable variants", y = "Count of individuals") +
  theme_minimal() +
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10, face = "bold")))

ggsave("actionable_variants.jpeg",
       plot = each_sample, device = "jpeg", height = 5, width = 6, units = "in", dpi = 600)

# Part B

(each_group <- ggplot(actionable, aes(x = as.factor(total_actionable), fill = Linguistic_tribal)) +
  geom_bar(color = "black") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  geom_text(stat = "count", aes(label = after_stat(count)), 
            vjust = -0.5, size = 3) +
  labs(x = "Number of actionable variants", 
       y = "Count of individuals") +
  facet_wrap(~ Linguistic_tribal, nrow = 1) +
  scale_fill_manual(values = c("#4682b4", "#8A8AF4", "#000080", "#FAC4C4", "#cd3333", "#C7F2C7", "#7ccd7c")) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10, face = "bold"),
    strip.text = element_text(size = 10, face = "bold", color = "white"),
    strip.background = element_rect(fill = "#4d4d4d", color = NA),
    legend.position = "none"
  ))

ggsave("actionable_variants_by_ling_tribal.jpeg",
       plot = each_group, device = "jpeg", height = 7, width = 15, units = "in", dpi = 600)
