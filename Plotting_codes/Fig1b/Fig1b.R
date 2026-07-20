library(tidyverse)
library(patchwork)

# Loading data ------------------------------------------------------------

#This input file is the supplementary table from the paper
source_data <- read.delim("./Supplementary_Table_S1.4.tsv", skip = 1)

df <- source_data %>%
  select(Numeric.Code, Number.of.samples, Tribe, Linguistic.Group) %>%
  mutate(Tribe = replace_values(Tribe,
                                "Yes" ~ "T",
                                "No"  ~ "NT",
                                 NA   ~ "CAO" )) %>%
  unite("group", Linguistic.Group, Tribe, remove = F) %>%
  mutate(group = replace_values(group, "NA_CAO" ~ "CAO")) %>%
  mutate(pop_bg_color = case_when(
    group == "IE_NT" ~ "#FAC4C4",
    group == "IE_T"  ~ "#CD3333",
    group == "DR_NT" ~ "#8A8AF4",
    group == "DR_T"  ~ "#000080",
    group == "AA_T"  ~ "#4682B4",
    group == "TB_NT" ~ "#C7F2C7",
    group == "TB_T"  ~ "#7CCD7C",
    group == "CAO"   ~ "#969696"
  )) %>%
  mutate(pop_font_color = case_when(
    group == "IE_T" ~ "#FFFFFF",
    group == "DR_T" ~ "#FFFFFF",
    group == "AA_T" ~ "#FFFFFF",
    .default = "#000000"
  )) %>%
  mutate(size_bg_color = rep(c("#F2F2F2", "#DEDEDE"), length.out = n()))

df$Numeric.Code <- factor(df$Numeric.Code, levels = df$Numeric.Code)


# Creating the plots ------------------------------------------------------

sample_sizes <- ggplot(df, aes(x = Numeric.Code, y = 1)) +
  geom_tile(aes(fill = size_bg_color)) +
  geom_text(aes(label = Number.of.samples), angle = 90, size = 9) +
  scale_fill_identity() +
  scale_x_discrete(expand = c(0, 0)) +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0))


population_ids <- ggplot(df, aes(x = Numeric.Code, y = 1, fill = pop_bg_color)) +
  geom_tile(width = 0.95, height = 0.1) +
  geom_text(aes(label = Numeric.Code, color = pop_font_color), size = 8, fontface = "bold") +
  scale_color_identity() +
  scale_fill_identity() +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0)) +
  theme(legend.position = "none")


# Putting the plots together and saving ----------------------

combined_plot <- sample_sizes / population_ids +
  plot_layout(heights = c(2.5, 1))

ggsave("Fig1b.png", width = 12000, height = 500, units = "px")
