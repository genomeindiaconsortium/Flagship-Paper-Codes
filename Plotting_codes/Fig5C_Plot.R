
data2=read.csv("MAF_Bins.txt",sep="\t")

#================This file looks like this===========================================

#MAF_Bins	      HRC_OVERLAP	  GASP_OVERLAP	  GI_OVERLAP	 TOPMED_OVERLAP
#[0,0.001)	    0.430299	    0.293192	      0.721493	   0.451242
#[0.001,0.002)	0.533775	    0.296517	      0.822272	   0.585903
#[0.002,0.003)	0.610961	    0.383502	      0.864347	   0.718665
#....

#<Panel>_OVERLAP : Mean panel-specific Rsq/INFO per MAF bin for overlapping variants.
#=====================================================================================

colors <- c("GI" = "forestgreen", "TOPMED" = "dodgerblue", "HRC" = "orange","GASP"="purple")

#Visualization
breaks_to_show <- c(
  "[0,0.001)", 
  "[0.009,0.01)",
  "[0.05,0.1)",
  "[0.45,0.5)"
)
upper_end <- sub(".*,(.*)\\)", "\\1", data2$MAF_Bins)
selected_bins <- which(data2$MAF_Bins %in% breaks_to_show)
ticks_df <- data.frame(
  x = 1:length(data2$MAF_Bins),                  # discrete positions
  bin = data2$MAF_Bins,
  label = ifelse(data2$MAF_Bins %in% breaks_to_show, upper_end, ""),
  tick_length = ifelse(data2$MAF_Bins %in% breaks_to_show, 0.05, 0.02),  # long vs short
  label_face = ifelse(data2$MAF_Bins %in% breaks_to_show, "bold", "plain")
)

#
gfg_plot4 = ggplot(data2, aes(x=MAF_Bins))+   
  geom_line(aes(y = GI_OVERLAP, color ="GI"), group=1,linewidth = 1.5)+
  geom_line(aes(y = TOPMED_OVERLAP, color = "TOPMED"), group=2,linewidth = 1.5) + 
  geom_line(aes(y = HRC_OVERLAP, color = "HRC"), group=3,linewidth = 1.5)+
  geom_line(aes(y = GASP_OVERLAP, color = "GASP"), group=4,linewidth = 1.5)+
  theme_minimal()+
  theme(panel.grid.minor = element_line(colour = "grey95"),panel.grid.major = element_line(colour = "grey95"))+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y = element_text(size=14, color="black"))+
  labs(x="Minor allele frequencies",y= "Imputation Accuracy Measure")+
  theme(legend.background = element_rect(fill = "white",colour="grey80"))+
  theme(legend.position = c(.89, .35),legend.justification = c("right", "bottom"),legend.box.just = "left",
        legend.title = element_blank(),legend.text = element_text(size=14),
        legend.key.size  = unit(1, "lines"),              # small key boxes
        legend.key.width = unit(1, "lines"),
        axis.title.y = element_text(size = 18),
        axis.title.x = element_text(size = 18),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())+
  scale_colour_discrete(labels = function(x) str_wrap(x, width = 5))+
  scale_color_manual(values=colors,breaks=c("GI","TOPMED","HRC","GASP"))+
  coord_cartesian(clip = "off") +
  geom_segment(data = ticks_df,
               aes(x = x, xend = x, y = -0.01, yend = -0.01 - tick_length),
               inherit.aes = FALSE,
               color = "black",
               size = 0.8) +  # tick thickness
  geom_text(data = ticks_df[ticks_df$label != "", ],
            aes(x = x, y = -0.09, label = label),
            inherit.aes = FALSE, hjust=0.5,
            size = 5)
gfg_plot4

ggsave("./R2vsMAF_Overlap_v1.png",
       gfg_plot4,
       width = 6,      # width in inches (or use units = "cm" etc)
       height = 6,     # height
       units = "in",   # units: "in", "cm" or "mm"
       dpi = 600       # resolution
)