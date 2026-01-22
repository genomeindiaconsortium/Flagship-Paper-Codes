data2=readxl::read_excel("Overlap_MAFvs_R2.xlsx",sheet=2)

#=======Input Files look like this==============================

#INFO	    MAF_Bins	    ALL_COUNT	  Panel	  MAF_range
#[0-0.2]	[0,0.001)	    387367	    TOPMED	Rare
#[0-0.2]	[0.009,0.01)	2812	      TOPMED	Rare
#...
#[0-0.2]	[0.01,0.015)	23456	      HRC	    Low Frequency
#...
#[0-0.2]	[0.05,0.1)	  36974	      GASP	  Common

#===============================================================

data2$INFO <- factor(data2$INFO, levels=c('[0-0.2]', '(0.2-0.4]', '(0.4-0.6]', '(0.6-0.8]','(0.8-1]'))
data2$MAF.range <- factor(data2$MAF_range,levels=c("Rare","LowFreq","Common"))
data2$MAF_group <- ifelse(data2$MAF_range=="Rare","MAF<1%","MAF>1%")

data2$INFO_new <- dplyr::case_when(
  data2$INFO %in% c('[0-0.2]', '(0.2-0.4]') ~ "0–0.4",
  data2$INFO == '(0.4-0.6]' ~ "0.4–0.6",
  data2$INFO == '(0.6-0.8]' ~ "0.6–0.8",
  data2$INFO == '(0.8-1]'   ~ "0.8–1",
  TRUE ~ NA_character_
)
data2$INFO_new <- factor(
  data2$INFO_new,
  levels = c("0–0.4","0.4–0.6","0.6–0.8","0.8–1")
)

data2$Panel <- ifelse(data2$Panel=="GASP","GAsP",data2$Panel)


data2$Panel=factor(data2$Panel, levels=c("GI","TOPMED","HRC","GAsP"))
data2$MAF_group <- as.factor(data2$MAF_group)

count2=ggplot(data2, aes(x = Panel, y = ALL_COUNT, fill = INFO_new)) +
  scale_fill_brewer(palette = "BrBG")+
  geom_col(position = position_stack())+
  facet_grid(~ MAF_group)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust=1,size=16),
        ,axis.text.y = element_text(size=16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size=16),legend.text = element_text(size=16),
        strip.text = element_text(size = 16, face = "bold"),
        plot.title = element_text(hjust=0.3, size=11))+
  labs(y="Counts",x="",fill = "Accuracy")+
  theme(
    legend.position = "right",       
    legend.direction = "vertical",   
    legend.justification = "left",
    legend.box.just = "left"
  ) +
  guides(fill = guide_legend(ncol = 1,byrow = TRUE )) 

ggsave("Count_Overlap_28112025.png",
       count2,
       width = 7,      
       height = 6,     
       units = "in",   
       dpi = 600       
)