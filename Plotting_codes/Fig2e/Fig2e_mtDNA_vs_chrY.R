library(dplyr)
library(ggplot2)
library(stringr)
library(ggrepel)



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
  labs(x = "Alpha Diversity of mtDNA Haplogroups", y = "Alpha Diversity of Y Haplogroups")

ggsave("./Fig2e.png", width = 12, height = 7, dpi = 300, bg = "white")