#FIGURE : 3D 
#ROH LENGTH ACROSS THE SAMPLED LINGUISTIC AND TRIBAL STATUS
#PACKAGE REQUIRED 
library(dplyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(tidyr)
library(ggh4x)  # this for the strip introduced in the plot



#reading the data

df <- read.csv("3D_ROH_withAJFINN.tsv", sep = "\t", header = T)
head(df)

metadata <- read.csv("3D_ROH_metadata_withAJFINN.tsv", sep = "\t", header = T)
head(metadata)

df <- df %>%
  mutate(ROH_bin = case_when(
    roh_len_cM < 10 ~ "1-10 cM",
    roh_len_cM <15 & roh_len_cM >= 10 ~ "10-15 cM",
    roh_len_cM >=15 ~ "≥15 cM"
  ))

roh_sum <- df %>%
  group_by(IID, ROH_bin) %>%
  summarise(sum_roh = sum(roh_len_cM), .groups = "drop")

roh_sum <- left_join(roh_sum, metadata, join_by(IID==IID))  
head(roh_sum)

roh_order <- roh_sum %>% group_by(LG_Tribe,IID) %>% 
  summarise(total_roh = sum(sum_roh), .groups="drop") %>%
  arrange(LG_Tribe,total_roh)


#plot asthetic settings :

roh_sum$IID <- factor(roh_sum$IID, levels = roh_order$IID)

roh_sum$ROH_bin = factor(roh_sum$ROH_bin, levels = c("≥15 cM", "10-15 cM", "1-10 cM"))

roh_sum$LG_Tribe <- roh_sum$LG_Tribe %>% replace_na("CAO")

roh_sum <- roh_sum %>%
  mutate(LG_Tribe_name = case_when(
    LG_Tribe == "AA_Tribe" ~ "AA_T",
    LG_Tribe == "DR_NonTribe" ~ "DR_NT",
    LG_Tribe == "DR_Tribe" ~ "DR_T",
    LG_Tribe == "IE_NonTribe" ~ "IE_NT",
    LG_Tribe == "IE_Tribe" ~ "IE_T",
    LG_Tribe == "TB_Tribe" ~ "TB_T",
    LG_Tribe == "TB_NonTribe" ~ "TB_NT",
    LG_Tribe == "AJ" ~ "AJ",
    LG_Tribe == "FIN" ~ "F\nI\nN",
    LG_Tribe ==  "CAO" ~ "C\nA\nO"
  ))



roh_sum$LG_Tribe_name <- gsub("_", "\n", roh_sum$LG_Tribe_name)
roh_sum$LG_Tribe_name <- factor(roh_sum$LG_Tribe_name, 
                                levels=c("F\nI\nN","AJ","AA\nT","DR\nNT","DR\nT","IE\nNT", "IE\nT","TB\nNT","TB\nT","C\nA\nO"))


#making the code :

# hex code acoording to what we have been using 
strip_colors <- c("FIN"=  "cyan", "AJ" = "#fff400","AA_T" ="#4682b4","DR_NT"="#8A8AF4","DR_T"="#000080" ,"IE_NT"="#FAC4C4", "IE_T"="#cd3333", "TB_T"="#7ccd7c","TB_NT"="#C7F2C7", "CAO" = "grey")

panel_widths <- roh_sum %>%
  count(LG_Tribe_name) %>%
  mutate(sqrt_n = sqrt(n)-5) %>%
  pull(sqrt_n)

ggplot(roh_sum, aes(x = IID, y = sum_roh, fill = ROH_bin)) +
  geom_bar(stat = "identity", width = 5) +
  labs(
    x = "Individuals",
    y = "Total ROH Length (cM)",
    fill = "ROH Bin"
  ) +
  scale_fill_manual(values = c("1-10 cM" ="#FEB24B",
                               "10-15 cM" = "#E3201C",
                               "≥15 cM"= "black")) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  
  facet_grid2(~ LG_Tribe_name, 
              scales = "free_x", 
              space = "free_x", 
              switch = "x",
              strip = strip_themed(
                background_x = elem_list_rect(fill = strip_colors)
              )
  ) + 
  force_panelsizes(cols = panel_widths) +
  
  theme_classic() +
  guides(
    fill = guide_legend(override.aes = list(shape = 21, size=8, color = "black"))
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 24),
    panel.border = element_blank(),
    plot.title = element_text(size = 24),
    axis.title = element_text(size = 24),
    strip.text.x = element_text(size = 0, lineheight = 1),
    legend.title = element_text(size = 20),
    legend.text = element_text(size=18),
    legend.position = "inside",
    legend.position.inside = c(0.6, 0.8),
    strip.placement = "outside"
  )


#SAVING THE PLOTS :
ggsave("./Fig3d.png", width=20, height=5, dpi=300)

