library(tidyverse)
library(sf)
library(ggimage)


# Loading data ------------------------------------------------------------

india <- read_sf("./India_map/polymap15m_area.shp") %>% st_transform(4326)

#This input file is the supplementary table from the paper
data <- read.delim("./Supplementary_Table_S1.4.tsv", skip = 1)  %>%
        mutate(Tribe = replace_values(Tribe, "Yes" ~ "T",
                                             "No"  ~ "NT",
                                              NA   ~ "CAO")) %>%
        unite("group", Linguistic.Group, Tribe, remove = F) %>%
        mutate(group = replace_values(group, "NA_CAO" ~ "CAO")) %>%
        mutate(point_shape = case_when(
          Tribe == "T"   ~ 22,
          Tribe == "NT"  ~ 21,
          Tribe == "CAO" ~ 24
        )) %>%
        mutate(point_color = case_when(
          group == "IE_NT" ~ "#FAC4C4",
          group == "IE_T"  ~ "#CD3333",
          group == "DR_NT" ~ "#8A8AF4",
          group == "DR_T"  ~ "#000080",
          group == "AA_T"  ~ "#4682B4", 
          group == "TB_NT" ~ "#C7F2C7",  
          group == "TB_T"  ~ "#7CCD7C",
          group == "CAO"   ~ "#969696"
        )) %>%
        mutate(font_color = case_when(
          group == "IE_T"  ~ "#FFFFFF",
          group == "DR_T"  ~ "#FFFFFF",
          group == "AA_T"  ~ "#FFFFFF",
          .default = "#000000"
        ))


# Plotting the map --------------------------------------------------------

gi <- ggplot(data) +
        geom_sf(data = india, fill = NA, color = "transparent") +
        geom_point(aes(x = Longitude, y = Latitude, shape = point_shape, 
                       fill = point_color), color = "black", size = 14, stroke = 2) +
        geom_text(aes(x = Longitude, y = Latitude, label = Numeric.Code, 
                      color = font_color), fontface = "bold", size = 7) +
        scale_fill_identity() +
        scale_shape_identity() +
        scale_color_identity() +
        theme_void() 


# Overlaying and saving ---------------------------------------------------

# The background map depicting linguistic group distributions, was manually created using external image editor
gi_bg <- ggbackground(gi, background = "background_map.png")

ggsave("Fig1a.png", gi_bg, width = 9530, height = 7760, units = "px")
