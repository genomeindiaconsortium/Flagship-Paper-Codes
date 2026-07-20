################################################################################
# Heatmap of Population-wise Pharmacogenomic Variant Frequencies
#
# This script generates the main heatmap figure from Supplementary Table S13.6.
#
# INPUT FILE:
# The input dataframe should be an Excel sheet (e.g., final_pgx_data.xlsx)
# derived from Supplementary Table S13.6.
#
# Each row should represent:
# One variant × one population combination.
#
# REQUIRED COLUMNS:
#
# 1. populations
#    - Population identifier.
#    - Should include:
#         1–82 for GenomeIndia populations
#         SAS, AFR, AMR, EAS, EUR for global populations
#
# 2. variant
#    - Pharmacogenomic variant identifier.
#    - Example: rs4244285
#
# 3. gene
#    - Gene corresponding to the variant.
#    - Example: CYP2C19
#
# 4. frequency
#    - Allele frequency of the variant in that population.
#    - Numeric column.
#
# 5. type
#    - Pharmacogenomic evidence category.
#    - Allowed values:
#         Dosage
#         Toxicity
#         Efficacy
#         Other
#
# 6. Testing
#    - CPIC/FDA testing recommendation category.
#    - Allowed values:
#         Testing Required
#         Testing Recommended
#         Actionable PGx
#
# 7. chemicals
#    - Drug or chemical name associated with the variant.
#    - Example: clopidogrel
#
# 8. category
#    - Therapeutic category of the drug.
#    - Example: Cardiovascular agents
#
# 9. Linguistic_group_Tribe
#    - Linguistic/tribal classification of each population.
#    - Allowed values:
#         AA_T
#         DR_T
#         DR_NT
#         IE_T
#         IE_NT
#         TB_T
#         NT
#         Global
#
# OUTPUT:
# The script produces a heatmap showing:
#   - Z-score normalized allele frequencies
#   - Population linguistic grouping
#   - Variant testing recommendation category
#   - Drug annotations
#
################################################################################

library(readxl)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(gridtext)

imp_var <- read_excel("final_pgx_data.xlsx")

imp_var <- imp_var |>
  mutate(populations = as.character(populations)) |>
  mutate(pop_num = suppressWarnings(as.numeric(populations))) |>
  arrange(is.na(pop_num), pop_num) |>
  mutate(populations = factor(populations, levels = c(as.character(1:82), "SAS", "AFR", "AMR", "EAS", "EUR")),
         type = factor(type, levels = c("Dosage", "Toxicity", "Efficacy", "Other")),
         Testing = factor(Testing, levels = c("Testing Required", "Testing Recommended", "Actionable PGx")),
         variant_gene = paste0(variant, " (<i>", gene, "</i>)"),
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

linguistic_colors <- c(
  "AA_T" = "#4682b4", "DR_T" = "#000080", "DR_NT" = "#8A8AF4",
  "IE_T" = "#cd3333", "IE_NT" = "#FAC4C4", "TB_T" = "#7ccd7c",
  "NT" = "#C7F2C7", "Global" = "#000000"
)
type_colors <- c("Toxicity"="#c6caca", "Dosage"="#c6caca", "Efficacy"="#c6caca", "Other"="#c6caca")

min_val <- min(imp_var$norm_frequency, na.rm = TRUE)
max_val <- max(imp_var$norm_frequency, na.rm = TRUE)

col_fun <- colorRamp2(
  c(min_val, -2.5, 0, 2.5, 5, max_val),
  c("#2E627AFF", "#1283c0ff", "#F5F5F5FF", "#D6604DFF", "#B2182BFF", "#67001FFF")
)

test_symbols <- c(
  "Testing Required" = "<span style='font-size:50pt;'>\u25CF</span>",
  "Testing Recommended" = "<span style='font-size:35pt;'>\u25B2</span>",
  "Actionable PGx" = "<span style='font-size:40pt;'>\u25A0</span>"
)

variant_labels <- paste0(test_symbols[row_meta$Testing], " ", rownames(mat))

row_meta <- imp_var %>%
  distinct(variant_gene, drug_category, type, Testing)
row_meta <- row_meta[match(rownames(mat), row_meta$variant_gene), ]

col_meta <- imp_var %>%
  distinct(populations, Linguistic_group_Tribe)
col_meta <- col_meta[match(colnames(mat), col_meta$populations), ]

col_ha <- HeatmapAnnotation(
  Linguistic_Group = col_meta$Linguistic_group_Tribe,
  col = list(Linguistic_Group = linguistic_colors),
  spacer = anno_empty(height = unit(1, "mm"), border = FALSE),
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
  
  Variant = anno_text(
    gt_render(variant_labels),
    gp = gpar(
      fontsize = 32,
      fontface = "bold",
      col = "black"
    ),
    just = "right",
    location = 0.99
  ),
  
  col = list(Type = type_colors),
  
  gap = unit(8, "mm"),
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
  row_title_gp = gpar(fontsize = 28, fontface = "bold"),
  column_split = col_meta$Linguistic_group_Tribe,
  column_title_gp = gpar(fontsize = 32, fontface = "bold"),
  column_names_gp = gpar(fontsize = 32, fontface = "bold"),
  row_gap = unit(4, "mm"),
  column_gap = unit(4, "mm"),
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
  title_gp = gpar(fontsize = 32, fontface = "bold"),
  labels_gp = gpar(fontsize = 28),
  legend_width = unit(16, "cm"),
  grid_height = unit(8, "mm"),
  title_gap = unit(5, "mm")
)

lg_testing <- Legend(
  labels = c(
    " Testing Required",
    " Testing Recommended",
    " Actionable PGx"
  ),

  title = "Testing",
  type = "points",
  pch = c(16, 17, 15),
  size = unit(8, "mm"),
  legend_gp = gpar(
    col = "black"
  ),
  title_gp = gpar(
    fontface = "bold",
    fontsize = 32
  ),
  labels_gp = gpar(
    fontsize = 32
  ),
  direction = "horizontal",
  title_gap = unit(5, "mm"),
  row_gap = unit(2, "mm")
)

combined_legends <- packLegend(
  lg_testing, lg_z, 
  direction = "horizontal",
  gap = unit(6, "cm")
)

png("final_heatmap.png", width = 60, height = 15, units = "in", res = 600)

draw(
  ht,
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom",
  merge_legend = FALSE,
  heatmap_legend_list = list(combined_legends)
)

dev.off()