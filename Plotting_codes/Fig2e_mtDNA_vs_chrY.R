library(dplyr)
library(ggplot2)
library(stringr)
library(ggrepel)



df_abdiversity <- read.table('./alpha_beta_admixture_diversity.txt', header = T, sep = "\t")
head(df_abdiversity)

colors <- c("#cd3333", "#000080", "#dbdb61", "#6b9bc3", "#96d796")
names(colors) <- c("ANI", "ASI", "ASI-D", "AAA", "ATB")

df_abdiversity %>%
  count(class)

df_abdiversity %>%
  mutate(ancestry = case_when(
    Population_Order %in% seq(1, 32) ~ "ANI", 
    Population_Order %in% seq(33, 60) ~ "ASI",
    Population_Order %in% seq(61, 62) ~ "ASI-D",
    Population_Order %in% seq(63, 74) ~ "AAA",
    Population_Order %in% seq(75, 82) ~ "ATB",
    Population_Order == "83" ~ "CAO",
    .default = "Others"
  )) %>%
  mutate(shape = case_when(
    Tribe_Shape == "Circle" ~ 21,
    Tribe_Shape == "Square" ~ 22,
    .default = 0
  )) %>%
  ggplot(aes(x = alpha_z, y = beta_z)) +
  geom_point(aes(fill = ancestry, shape = shape), stroke = 1, size = 5) +
  theme_classic(24) +
  theme(axis.text = element_text(color = "black")) +
  geom_text_repel(aes(label = Population_Order), point.padding = 1, box.padding = 0.4, max.overlaps = 50, min.segment.length = 0.5, size = 5) + 
  scale_fill_manual(values = colors, guide = F) +
  scale_shape_identity() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  annotate("rect", xmin = -3, xmax = -2, ymin = -1.8, ymax = -0.1, fill = "steelblue", alpha = 0.15) +
  annotate("rect", xmin = -1.9, xmax = -0.45, ymin = -1.5, ymax = 0.5, fill = "brown1", alpha = 0.1) +
  annotate("rect", xmin = 1.35, xmax = 1.8, ymin = -1, ymax = 1.5, fill = "darkgoldenrod1", alpha = 0.2) +
  annotate("rect", xmin = -0.4, xmax = 0.7, ymin = 1.2, ymax = 2, fill = "darkslategray4", alpha = 0.2) +
  annotate("rect", xmin = -1, xmax = 0.6, ymin = 2.7, ymax = 3.5, fill = "darkgray", alpha = 0.2) +
  labs(x = "Alpha Diversity of Estimated Ancestry", y = "Beta Diversity of Estimated Ancestry")

ggsave("./Fig2d.png", width = 12, height = 7, dpi = 300, bg = "white")  



df <- read.table("./mtDNA_vs_chrY_AlphaBetaDiversity.tsv", sep = "\t", header = T, stringsAsFactors = F)

head(df)


colors <- c("#cd3333", "#000080", "#dbdb61", "#6b9bc3", "#96d796")
names(colors) <- c("ANI", "ASI", "ASI-D", "AAA", "ATB")

df %>%
  mutate(ancestry = case_when(
    Population.Order %in% seq(1, 32) ~ "ANI", 
    Population.Order %in% seq(33, 60) ~ "ASI",
    Population.Order %in% seq(61, 62) ~ "ASI-D",
    Population.Order %in% seq(63, 74) ~ "AAA",
    Population.Order %in% seq(75, 82) ~ "ATB",
    Population.Order == "83" ~ "CAO",
    .default = "Others"
  )) %>%
  mutate(shape = if_else(str_detect(Ethnolinguistic.Group, "NonTribe"), 21, 22)) %>%
  filter(ancestry != "CAO") %>%
  ggplot(aes(x = mtDNA.Diversity..Z.score., y = Y.Diversity..Z.score.)) +
  geom_point(aes(fill = ancestry, shape = shape), stroke = 1, size = 5) +
  geom_text_repel(aes(label = Population.Order), point.padding = 1, box.padding = 0.4, max.overlaps = 50, min.segment.length = 0.5, size = 5) + 
  theme_classic(24) +
  theme(axis.text = element_text(color = "black")) +
  scale_fill_manual(values = colors, guide = F) +
  scale_shape_identity() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  annotate("rect", xmin = -2.5, xmax = -0.3, ymin = -3.1, ymax = -1.1, fill = "coral4", alpha = 0.1) +
  annotate("rect", xmin = 1.1, xmax = 2, ymin = 0.3, ymax = 1, fill = "darkgreen", alpha = 0.15) +
  annotate("rect", xmin = 1, xmax = 1.7, ymin = -2.5, ymax = -1.5, fill = "cyan4", alpha = 0.1) +
  xlim(-2.6, 2.5) + ylim(-3.1, 2) +
  labs(x = "Alpha Diversity of mtDNA haplogroups", y = "Alpha Diversity of Y haplogroups")

ggsave("./Fig2e.png", width = 12, height = 7, dpi = 300, bg = "white")
