
rm(list = ls())
library(ggplot2)

data <- data.frame(
  AF_bin = factor(rep(c("Singleton", "Doubletons", "<0.1", ">=0.1 to <1", ">=1 to <5", ">=5"), each = 4), 
                  levels = c("Singleton", "Doubletons", "<0.1", ">=0.1 to <1", ">=1 to <5", ">=5")),
  Variant_Type = rep(c("SNV_Total", "SNV_Novel", "INDEL_Total", "INDEL_Novel"), times = 6),
  Value = c(
    55951809, 29986045, 3240269, 2107007, # Singleton
    16562454, 6022615, 1033605, 474177,   # Doubletons
    32155953, 4770340, 2548505, 520387,   # MAC3-UltraRare
    8836252, 19377, 732320, 27546,        # MAC3-Rare
    2540788, 12, 213867, 18388,           # MAC3-LowFreq
    5692113, 0, 430954, 88447            # MAC3-Common
  )
)

# Add separate columns for "Variant" and "Status"
data$Variant <- gsub("_.*", "", data$Variant_Type)   # Extract SNV or INDEL
data$Status <- gsub(".*_", "", data$Variant_Type)    # Extract Total or Novel

# Filter data for SNV and INDEL separately
data_snv <- subset(data, Variant == "SNV")
data_indel <- subset(data, Variant == "INDEL")

max_y <- max(data$Value)
# Add the desired font
font_add_google("Lora", "safira") 
showtext_auto()

# Create SNV plot
plot_snv <- ggplot(data_snv, aes(x = AF_bin, y = Value, fill = Status)) +
  geom_bar(stat = "identity", position = position_dodge(width = 1), color = "black") +
  geom_text(aes(label = paste0(round(Value / 1e6, 2))), 
            position = position_dodge(width = 0.8), size = 6,vjust = -0.3, hjust = 0.5) +
  scale_fill_manual(values = c("Total" = "#1f78b4", "Novel" = "#a6cee3"), name = "SNV") +
  scale_y_continuous(limits = c(0, max_y),labels = function(x) paste0(x / 1e6), 
                     expand = expansion(mult = c(0, 0.6))) +  
  labs(   
    x = "MAF",
    y = "Number of Variants \n (million)"
  ) +
  theme_minimal() +
  theme(
    text = element_text(family = "safira"),  # Apply custom font
    axis.text.x = element_text(size = 20, colour = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 24, colour = "black"),
    axis.title.x = element_text(size = 24, face = "bold"),
    axis.title.y = element_text(size = 24, face = "bold"),
    legend.title = element_text(size = 24, face = "bold"),
    legend.text = element_text(size = 20, colour = "black"),
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
    legend.position = c(0.8, 0.8)
  )+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.justification = "center")

# Create INDEL plot
plot_indel <- ggplot(data_indel, aes(x = AF_bin, y = Value, fill = Status)) +
  geom_bar(stat = "identity", position = position_dodge(width = 1), color = "black") +
  geom_text(aes(label = paste0(round(Value / 1e6, 2))), 
            position = position_dodge(width = 0.8), vjust = -0.3, hjust = 0.5, size = 6) +
  scale_fill_manual(values = c("Total" = "#33a02c", "Novel" = "#b2df8a"), name = "INDEL") +
  scale_y_continuous(labels = function(x) paste0(x / 1e6), expand = expansion(mult = c(0, 0.5))) +  

  labs(
    x = "MAF bins",
    y = "Number of Variants \n (million)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 20, colour = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 24, colour = "black"),
    axis.title.x = element_text(size = 24, face = "bold"),
    axis.title.y = element_text(size = 24, face = "bold"),
    legend.title = element_text(size = 24, face = "bold"),
    legend.text = element_text(size = 20, colour = "black"),
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
    legend.position = c(0.8, 0.8)
  )+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.justification = "center")


spacer <- ggplot() + 
  theme_void() + 
  theme(plot.margin = margin(t = 50, b = 0, l = 0, r = 0)) 

combined_plot <- plot_grid(
  spacer,                                          
  plot_grid(plot_snv, plot_indel, ncol = 2,       
            align = "h", rel_widths = c(1.1, 1.1)),
  ncol = 1, rel_heights = c(0.1, 1)                
)


tiff("Figure2B.tiff", units="in", width=18, height=6, res=300)

print(combined_plot)

dev.off()
