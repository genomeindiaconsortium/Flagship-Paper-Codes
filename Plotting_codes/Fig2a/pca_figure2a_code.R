#PCA PLOT 
#FIGURE 2a

#REQUIRED PACKAGES 

library(ggplot2)
library(ggpubr)

#READING THE DATA :

pl_pca_fin <- read.csv("pca_2a/pca_2a_data.csv" , header = T )


#MAKING THE PLOT 
p4 <-ggplot(pl_pca_fin, aes(x = PC1, y = -(PC2), color = new_col , shape = Tribe) ) +
  geom_point(size = 1) + 
  scale_shape_manual(values = c(
    "Yes" = 15,
    "No" =16
  )) +  
  scale_color_manual(values = c(
    "IE_Yes" = "#cd3333",
    "DR_Yes" = "#000080",
    "AA_Yes" = "#4682b4",
    "TB_Yes" = "#7ccd7c",
    "IE_No" = "#FAC4C4",
    "DR_No" = "#8A8AF4",
    "AA_No" = "#99C5EC",
    "TB_No" = "#C7F2C7"
  )) +
  labs(
    x = "PC1",
    y = "PC2",
  ) +
  theme_pubr() +  
  theme( legend.position = "none")  # Remove legend )

print(p4)




