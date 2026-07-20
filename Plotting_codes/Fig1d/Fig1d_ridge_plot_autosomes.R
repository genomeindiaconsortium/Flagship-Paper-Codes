library(ggridges)
library(ggplot2)
library(dplyr)
library(patchwork)
library(ggpubr)
library(tidyverse)
library(showtext)


showtext_auto()
GI_data <- read.delim("/home/Sreelekshmi/Downloads/dec_2024_vcf_stats_all_chr_2.tsv")

GI_data_all <- GI_data
GI_data_all$Groups <- "All"

my_theme <- theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
                  axis.text = element_text(size = 16, colour = "black"),
                  axis.title = element_text(size = 20, color = "black"),
                  panel.grid.minor = element_blank(),
                  legend.text = element_text(size = 20, colour = "black"),
                  legend.title = element_blank(),
                  #legend.title = element_text(size = 20, colour = "black"),
                  # legend.key.size = unit(0.8, 'cm'),
                  axis.line = element_line(colour = "black")
                  )

GI_data_combined <- rbind(GI_data, GI_data_all) %>%
                    mutate (Groups = recode(Groups, "TB_2"="TB_NT", "TB_1"="TB_T", "DR_1"="DR_T", "AA_1"="AA_T", 
                                            "IE_1"="IE_T", "DR_2"="DR_NT", "IE_2"="IE_NT"))

GI_data_combined$Groups <- factor(GI_data_combined$Groups, levels = rev(c("All","AA_T","DR_NT","DR_T","IE_NT","IE_T","TB_NT","TB_T","CAO")))

#group_colors <- c("All" = "orange", "AA_T" = "#048BA8", "CAO" = "#FF579F", "DR_NT" = "#9DD9D2", "DR_T" = "#FFE900",
#                  "IE_NT" = "#4E148C", "IE_T"= "#904E55", "TB_NT"= "#DFBBB1", "TB_T" = "#119822")

group_colors <- c("All" = "orange", "AA_T" = "#4682B4", "CAO" = "grey45", "DR_NT" = "#8A8AF4", "DR_T" = "#1d1d89" ,
                  "IE_NT" = "#FAC4C4", "IE_T"= "#cf3b3b", "TB_NT"= "#C7F2C7", "TB_T" = "#7ccd7c")


median_all_SNVs <- GI_data_combined %>% filter(Groups == "All") %>% summarise(med = median(SNVs, na.rm = TRUE)) %>% pull(med)

plot2 <- ggplot(GI_data_combined, aes(y = Groups, x = SNVs, fill = Groups)) +
  geom_density_ridges(scale=1.3, quantile_lines = F, quantile_fun=function(SNVs,...)median(SNVs)) +
  geom_segment(aes(x = median_all_SNVs, xend = median_all_SNVs, y = 1, yend = 10.5),
               linetype = "dashed", size = 0.6, colour = "black") +
  theme_minimal() +
  labs(x = "SNVs", y = "Groups") +
  my_theme + 
  scale_y_discrete(expand = expansion(add = c(.3, 1.6))) +
  scale_x_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale())) + 
  scale_fill_manual(values = group_colors) + guides(fill = guide_legend(nrow = 2, byrow = FALSE, reverse = TRUE)) 
plot2

median_all_INDELs <- GI_data_combined %>% filter(Groups == "All") %>% summarise(med = median(INDELs, na.rm = TRUE)) %>% pull(med)
plot3 <- ggplot(GI_data_combined, aes(y = Groups, x = INDELs, fill = Groups)) +
  geom_density_ridges(scale=1.3, quantile_lines = F, quantile_fun=function(INDELs,...)median(INDELs)) +
  geom_segment(aes(x = median_all_INDELs, xend = median_all_INDELs, y = 1, yend = 10.5),
               linetype = "dashed", size = 0.6, colour = "black") +
  theme_minimal() +
  labs(x = "INDELs", y = "Groups") + 
  my_theme + 
  scale_y_discrete(expand = expansion(add = c(.3, 1.6))) +
  scale_x_continuous(breaks = c(220000, 240000, 260000, 280000),
                     labels = scales::label_number(scale_cut = scales::cut_short_scale())) + 
  scale_fill_manual(values = group_colors) +  guides(fill = guide_legend(nrow = 2, byrow = FALSE, reverse = TRUE))
plot3


median_all_Singletons <- GI_data_combined %>% filter(Groups == "All") %>% summarise(med = median(Singletons, na.rm = TRUE)) %>% pull(med)
plot4 <- ggplot(GI_data_combined, aes(y = Groups, x = Singletons, fill = Groups)) +
  geom_density_ridges(scale=1.3, quantile_lines = F, quantile_fun=function(Singletons,...)median(Singletons)) +
  geom_segment(aes(x = median_all_Singletons, xend = median_all_Singletons, y = 1, yend = 10.5),
               linetype = "dashed", size = 0.6, colour = "black") +
  theme_minimal() +
  labs(x = "Singletons", y = "Groups") +
  my_theme + 
  scale_y_discrete(expand = expansion(add = c(.3, 1.6))) +
  scale_x_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale())) + 
  scale_fill_manual(values = group_colors) + guides(fill = guide_legend(nrow = 2, byrow = FALSE, reverse = TRUE))

plot4

plot_grid <- plot2 + plot3 + plot4
plot_grid + plot_layout(guides = "collect", axes = "collect_y", nrow = 1) &  theme(legend.position = "bottom")

showtext_opts(dpi = 300)
ggsave(file=paste0("/home/Sreelekshmi/Downloads/pop_ridge_autosomes_new_6.png"), width = 12, height = 5.5, bg="white", dpi = 300)
showtext_opts(dpi = 96)
showtext_auto()


#############################################################################################################################

GI_data_legend <- GI_data %>% mutate (Groups = recode(Groups, "TB_2"="TB_NT", "TB_1"="TB_T", "DR_1"="DR_T", "AA_1"="AA_T", 
                          "IE_1"="IE_T", "DR_2"="DR_NT", "IE_2"="IE_NT"))

GI_data_legend$Groups <- factor(GI_data_legend$Groups, levels = rev(c("DR_NT","DR_T","IE_NT","IE_T","TB_NT","TB_T","AA_T","CAO")))


plot5 <- ggplot(GI_data_legend, aes(y = Groups, x = SNVs, fill = Groups)) +
  geom_density_ridges(quantile_lines = F, quantile_fun=function(SNVs,...)median(SNVs)) +
  geom_segment(aes(x = median_all_SNVs, xend = median_all_SNVs, y = 1, yend = 10.5),
               linetype = "dashed", size = 0.6, colour = "black") +
  theme_minimal() +
  labs(x = "SNVs", y = "Groups") +
  my_theme +  theme(legend.key.spacing.x = unit(0.6, "cm"), legend.key.spacing.y = unit(0.2, "cm")) +
  scale_y_discrete(expand = expansion(add = c(.3, 1.6))) +
  scale_x_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale())) + 
  scale_fill_manual(values = group_colors) + guides(fill = guide_legend(nrow = 2, byrow = FALSE, reverse = TRUE)) 
plot5

legend <- get_legend(plot5)
as_ggplot(legend)

showtext_opts(dpi = 300)
ggsave(file=paste0("/home/Sreelekshmi/Downloads/pop_ridge_autosomes_legend_2.png"), width = 6, height = 2, bg="white", dpi = 300)
showtext_opts(dpi = 96)
showtext_auto()
