library(dplyr)
library(ggplot2)
library(stringr)
library(ggrepel)
library(forcats)
library(paletteer)
library(patchwork)

prs_df <- read.table("./prs_input.tsv", header = T, sep = "\t")

head(prs_df)

paletteer_d("ggsci::nrc_npg")

colors <- c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF")
colors <- setNames(colors, c("GI", "UKBB_African", "UKBB_British", "UKBB_Indian"))

trait_order <- c("Height", "Weight", "BMI")
pop_order <- c("UKBB_British", "UKBB_Indian", "GI", "UKBB_African")

r2_plot <- prs_df %>%
  mutate(pop = if_else(Cohort == "GI", Cohort, paste0(Cohort, "_", Population))) %>%
  ggplot(aes(x = r2_med, y = fct_relevel(Trait, rev(trait_order)), color = fct_relevel(pop, rev(pop_order)))) +
  geom_point(position = position_dodge(width = 0.75), size = 4) +
  geom_hline(yintercept = c(1.5, 2.5), color = "black") +
  geom_errorbar(aes(xmin = r2_min, xmax = r2_max), position = position_dodge(width = 0.75), width = 0.4, linewidth = 1) +
  # geom_vline(xintercept = 0, linetype = "dashed") +
  theme_classic(20) +
  theme(axis.text = element_text(color = "black")) +
  labs(x = "Variance Explained", y = "Trait", color = "Cohort") +
  scale_color_manual(values = colors, guide = guide_legend(reverse = T))
  
r2_plot

ci_plot <- prs_df %>%
  mutate(pop = if_else(Cohort == "GI", Cohort, paste0(Cohort, "_", Population))) %>%
  ggplot(aes(x = coeff_median, y = fct_relevel(Trait, rev(trait_order)), color = fct_relevel(pop, rev(pop_order)))) +
  geom_point(position = position_dodge(width = 0.75), size = 4) +
  geom_hline(yintercept = c(1.5, 2.5), color = "black") +
  geom_errorbar(aes(xmin = coeff_min, xmax = coeff_max), position = position_dodge(width = 0.75), width = 0.4, linewidth = 1) +
  # geom_vline(xintercept = 0, linetype = "dashed") +
  xlim(0, 1.7) +
  theme_classic(20) +
  theme(axis.text = element_text(color = "black")) +
  labs(x = "Effect Size", y = "Trait", color = "Cohort") +
  scale_color_manual(values = colors, guide = guide_legend(reverse = T))
ci_plot

r2_plot + ci_plot + plot_layout(guides = "collect", axes = "collect_y")

r2_plot + ci_plot + plot_layout(guides = "collect", axes = "collect_y") &
  theme(legend.position = "bottom", legend.margin = margin(r = 2, unit = "cm"))

showtext_opts(dpi = 300)
ggsave("./PRS_combined.png", width = 9, height = 6.2, dpi = 300, bg="white")
showtext_opts(dpi = 96)
showtext_auto()
