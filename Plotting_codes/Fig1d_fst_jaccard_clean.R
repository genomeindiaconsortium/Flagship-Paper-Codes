library(tidyverse)
library(RColorBrewer)
library(showtext)
library(ggtext)

data <- read.csv("./fst_jacard_data.csv") %>%
  filter(Pop1 != "CAO" & Pop2 != "CAO") %>%
  rowwise() %>%
  mutate(
    PopA = min(Pop1, Pop2), 
    PopB = max(Pop1, Pop2)
  ) %>%
  ungroup() %>%
  distinct(PopA, PopB, .keep_all = TRUE)

all_pops <- sort(unique(c(data$PopA, data$PopB)))

# Common theme function
common_theme <- function() {
  theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      plot.margin = margin(10, 10, 10, 10),
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.background = element_rect(fill = "transparent", color = NA),
      legend.background = element_rect(fill = "transparent", color = NA),
      legend.box.background = element_rect(fill = "transparent", color = NA),
      legend.key = element_rect(fill = "transparent", color = NA)
    )
}

# FST Plot
plot_data_fst <- data %>%
  filter(PopA <= PopB) %>%
  mutate(
    PopA = factor(PopA, levels = rev(all_pops)),
    PopB = factor(PopB, levels = all_pops)
  )

fst_plot <- ggplot(plot_data_fst, aes(x = PopB, y = PopA, fill = NIBMG.EIGENSOFT.)) +
  geom_tile(color = "white", linewidth = 0.3, na.rm = TRUE) +
  geom_text(
    aes(
      label = ifelse(is.na(NIBMG.EIGENSOFT.), "", sprintf("%.3f", NIBMG.EIGENSOFT.)),
      color = NIBMG.EIGENSOFT. > 0.02
    ),
    size = 8,
    na.rm = TRUE
  ) +
  scale_fill_gradientn(
    colors = brewer.pal(9, "YlOrRd"),
    na.value = "white",
    limits = c(0, max(plot_data_fst$NIBMG.EIGENSOFT., na.rm = TRUE))
  ) +
  scale_color_manual(
    values = c("TRUE" = "white", "FALSE" = "black"),
    guide = "none"
  ) +
  labs(x = NULL, y = NULL, title = NULL, fill = expression("Weighted F"["ST"])) +
  common_theme() +
  theme(
    axis.text.x.top = element_text(
      angle = 45, hjust = 0, vjust = 0, 
      size = 22, color = "black", face = "plain"
    ),
    axis.text.y.right = element_text(
      hjust = 0, margin = margin(l = 5), 
      size = 22, color = "black", face = "plain"
    ),
    legend.position = "right",
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 20)
  ) +
  scale_x_discrete(position = "top") +
  scale_y_discrete(position = "right") +
  coord_fixed()

fst_plot

showtext_opts(dpi = 300)
ggsave("./FST.png", plot = fst_plot, width = 10, height = 8, dpi = 300, bg = "transparent")
showtext_opts(dpi = 96)
showtext_auto()

# Jaccard Plot
plot_data_jaccard <- data %>%
  mutate(
    PopA = factor(PopA, levels = all_pops),
    PopB = factor(PopB, levels = rev(all_pops))
  )

jaccard_plot <- ggplot(plot_data_jaccard, aes(x = PopA, y = PopB, fill = jacard_1)) +
  geom_tile(
    aes(fill = ifelse(jacard_1 == 0, NA, jacard_1)),
    color = "white", linewidth = 0.3, na.rm = TRUE
  ) +
  geom_text(
    aes(
      label = ifelse(jacard_1 != 0, sprintf("%.3f", jacard_1), NA_character_),
      color = jacard_1 >= 0.480
    ),
    size = 8,
    na.rm = TRUE
  ) +
  scale_fill_gradientn(
    name = "Jaccard index",
    colors = brewer.pal(9, "YlGnBu"),
    na.value = "white",
    limits = c(0.380, 0.600),
    labels = function(x) sprintf("%.3f", x)
  ) +
  scale_color_manual(
    values = c("TRUE" = "white", "FALSE" = "black"),
    guide = "none"
  ) +
  labs(x = NULL, y = NULL, title = NULL) +
  common_theme() +
  theme(
    axis.text.x.bottom = element_text(
      angle = 45, hjust = 0, vjust = 0, 
      size = 22, color = "black"
    ),
    axis.text.y.left = element_text(
      hjust = 0, margin = margin(l = 5), 
      size = 22, color = "black"
    ),
    legend.position = "left",
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 20),
    legend.key.height = unit(1, "cm")
  ) +
  scale_x_discrete(position = "bottom") +
  scale_y_discrete(position = "left") +
  coord_fixed()

jaccard_plot

showtext_opts(dpi = 300)
ggsave("./jacard.png", plot = jaccard_plot, width = 10, height = 8, dpi = 300, bg = "transparent")
showtext_opts(dpi = 96)
showtext_auto()
