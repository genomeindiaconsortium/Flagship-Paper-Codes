#PCA PLOT
#figure 2b 

#REQUIRED PACKAGES:

library(ggplot2)
library(ggpubr)
library(ggrepel)
library(dplyr)


#read the data :

mean_df2 <- read.csv("pca_2b/pca_2b_data.txt")


p <- ggplot(mean_df2, aes(x = pc1, y = -pc2)) +
  geom_point(aes(shape = Tribe,
                 fill = Linguistic_Group),
             size = 2.5,
             alpha = 0.8,
             color = "black",   # outline color
             stroke = 1) +      # outline thickness
  scale_shape_manual(values = c("GI_Yes" = 22,   
                                "GI_No" = 21,    
                                "GAsP"  = 24,    
                                "HGDP"  = 23)) + 
  scale_fill_manual(values = c(
    "IE_T"="#cd3333",
    "DR_T"="#000080",
    "AA_T" ="#4682b4",
    "TB_T"="#7ccd7c",
    "IE_NT"="#FAC4C4",
    "DR_NT"="#8A8AF4",
    "AA_NT" ="#99C5EC",
    "TB_NT"="#C7F2C7",
    "NEA"               = "orange",
    "SEA"               = "slategray1",
    "EAST_ASIA"         = "darkgrey",
    "CSA"= "purple"
  )) +
  labs(x = "PC1",
       y = "PC2",
       fill = "Linguistic Group / Region",
       shape = "Consortium Name / Tribe Status") +
  theme_pubr()+ theme(legend.position = "none")

print(p)






