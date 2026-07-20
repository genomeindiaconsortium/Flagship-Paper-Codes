# Load libraries
library(readxl)
library(ggplot2)
library(ggrepel)

# Read data
df <- read_excel(
  "PATH/Supplementary Tables.xlsx",
  sheet = "Table S7.1",
  skip = 2
)

# Define color palette
colorset <- c(
  "IE_T"  = "brown3",
  "DR_T"  = "navy",
  "AA_T"  = "steelblue",
  "TB_T"  = "palegreen3",
  "IE_NT" = "#FAC4C4",
  "DR_NT" = "#8A8AF4",
  "AA_NT" = "#99C5EC",
  "TB_NT" = "#C7F2C7",
  "CAO"   = "grey70",
  "AJ"    = "#fff400",
  "FIN"   = "cyan"
)


shapeset <- c(
  "NonTribe" = 21,
  "Tribe"    = 22,
  "CAO"      = 24
)

# Set factor order
df$`Ethnolinguistic Group` <- factor(
  df$`Ethnolinguistic Group`,
  levels = c(
    "IE_T", "IE_NT",
    "DR_T", "DR_NT",
    "AA_T", "AA_NT",
    "TB_T", "TB_NT",
    "CAO", "AJ", "FIN"
  )
)

# Create shape variable
df$shape_var <- ifelse(
  df$`Ethnolinguistic Group` == "CAO",
  "CAO",
  ifelse(df$Tribe == "Yes", "Tribe", "NonTribe")
)

# Plot
png("Figure3C.png", width = 4000, height = 3000, res = 400)

ggplot(df) +
  
  # Main points
  geom_point(
    aes(
      x = `median of LROH (cM)`,
      y = `median of NROH`,
      fill = `Ethnolinguistic Group`,
      shape = shape_var
    ),
    size = 4,
    color = "black",
    stroke = 0.7
  ) +
  
  # Highlight AJ and FIN
  geom_point(
    data = subset(df,
                  `Ethnolinguistic Group` %in% c("AJ", "FIN")),
    aes(
      x = `median of LROH (cM)`,
      y = `median of NROH`,
      fill = `Ethnolinguistic Group`,
      shape = shape_var
    ),
    size = 6,
    color = "black",
    stroke = 1
  ) +
  
  # Labels for all populations
  geom_text_repel(
    aes(
      x = `median of LROH (cM)`,
      y = `median of NROH`,
      label = `Numeric Code`
    ),
    size = 2.5,
    max.overlaps = 100
  ) +
  

  
  # Manual scales
  scale_fill_manual(values = colorset) +
  scale_shape_manual(values = shapeset) +
  
  # Axis labels
  labs(
    x = "Average Length of ROH (cM)",
    y = "Number of ROH"
  ) +
  
  # Theme
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.title = element_text(size = 20, face = "bold"),
    axis.text = element_text(size = 16),
    plot.title = element_text(face = "bold")
  )

dev.off()

