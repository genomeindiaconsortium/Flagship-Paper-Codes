library(readxl)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)


imp_var <- read.table("./final_pgx_data.txt", header = T, sep = "\t")

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

bio_geo_data <- imp_var %>% distinct(populations, Biogeography)
bio_geo_data <- bio_geo_data[match(colnames(mat), bio_geo_data$populations), ]

row_meta <- imp_var %>%
  distinct(variant_gene, drug_category, type, Testing)
row_meta <- row_meta[match(rownames(mat), row_meta$variant_gene), ]

col_meta <- imp_var %>%
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
                   gp = gpar(fontsize = 30, fontface = "italic"),
                   just = "left",
                   location = 0.01)
)

left_ha <- rowAnnotation(
  Type = row_meta$type,
  Variant = anno_text(rownames(mat), 
                      gp = gpar(fontsize = 28, 
                                col = testing_colors[row_meta$Testing]),
                      just = "right",
                      location = 0.99),
  col = list(Type = type_colors),
  gap = unit(6, "mm"), 
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
  # row_names_side = "left",
  row_split = row_meta$type,
  row_title_gp = gpar(fontsize = 20),
  column_split = col_meta$Linguistic_group_Tribe,
  column_title_gp = gpar(fontsize = 36),
  column_names_gp = gpar(
    fontsize = 30,
    col = bio_geo_colors[col_meta$Biogeography]
  ),
  # row_names_gp = gpar(...),
  row_gap = unit(5, "mm"),
  column_gap = unit(5, "mm"),
  top_annotation = col_ha,
  right_annotation = right_ha,
  left_annotation = left_ha,
  border = TRUE,
  rect_gp = gpar(col = "white", lwd = 1),
  heatmap_legend_param = list(
    title = "Z-score Frequency",
    at = seq(floor(min_val), ceiling(max_val), by = 2),
    labels = seq(floor(min_val), ceiling(max_val), by = 2),
    title_gp = gpar(fontsize = 36),
    labels_gp = gpar(fontsize = 28),
    direction = "horizontal"
  )
)

lg_testing <- Legend(
  direction = "horizontal",
  labels = names(testing_colors),
  title = "Testing",
  legend_gp = gpar(fill = testing_colors),
  title_gp = gpar(fontsize = 36),
  labels_gp = gpar(fontsize = 30),
)

showtext_opts(dpi = 600)
png("./final_heatmap.png",
    width = 50, height = 15, units = "in", res = 600)
draw(ht, 
     heatmap_legend_side = "bottom",
     annotation_legend_side = "bottom",
     merge_legend = TRUE,
     annotation_legend_list = list(lg_testing))
dev.off()
showtext_opts(dpi = 96)
showtext_auto()
