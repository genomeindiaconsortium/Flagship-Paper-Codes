library(ggplot2)
library(extrafont)
library(showtext)
font_add_google("PT Sans Narrow", "pt_sans_narrow")

showtext_auto()

GI_barplot_Fig1_SNV <- read.delim("./Fig1_novel_var_barplot_SNV.tsv")

GI_barplot_Fig1_SNV$MAF <- factor(GI_barplot_Fig1_SNV$MAF, levels = c("AC=1","Ultra Rare","Rare","Low","Common"))
GI_barplot_Fig1_SNV$SNV <- factor(GI_barplot_Fig1_SNV$SNV, levels = c("Total","Novel"))

my_theme <- theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
                  axis.text = element_text(size = 16, colour = "black"),
                  axis.title = element_text(size = 20, color = "black"),
                  panel.grid.minor = element_blank(),
                  legend.text = element_text(size = 20, colour = "black"),
                  legend.title = element_text(size = 20, colour = "black"),
                  legend.position = c(0.8,0.8),
                  axis.line = element_line(colour = "black")
                  )
                                        
p1 <- ggplot(GI_barplot_Fig1_SNV, aes(x=MAF, y=Count, fill=SNV)) + 
      geom_bar(width=0.7,position=position_dodge(width=0.8), stat="identity",colour="black") +
      theme_classic() +
      labs(x = "MAF", y = "Number of variants\n(millions)") +
      my_theme + geom_text(aes(label = ifelse(Count < 1, Count, sprintf("%.1f", Count))), family = "pt_sans_narrow",  
                           position = position_dodge(width = 0.8), color = "black", size = 6, vjust = -0.3) +
      scale_y_continuous(limits = c(0, 65), breaks = seq(0, 65, 10), expand = c(0, 0)) +
      scale_fill_manual(values = c("Total"="#c51b7d","Novel"="#f1b6da"))

#"#c51b7d", "#f1b6da"

#"#762a83", "#c2a5cf"

p1 

GI_barplot_Fig1_INDEL <- read.delim("./Fig1/Fig1_novel_var_barplot_INDEL.tsv")

GI_barplot_Fig1_INDEL$MAF <- factor(GI_barplot_Fig1_INDEL$MAF, levels = c("AC=1","Ultra Rare","Rare","Low","Common"))
GI_barplot_Fig1_INDEL$INDEL <- factor(GI_barplot_Fig1_INDEL$INDEL, levels = c("Total","Novel"))

p2 <- ggplot(GI_barplot_Fig1_INDEL, aes(x=MAF, y=Count, fill=INDEL)) + 
  geom_bar(width=0.7,position=position_dodge(width=0.8), stat="identity",colour="black") +
  theme_classic() +
  labs(x = "MAF", y = "Number of variants\n(millions)") +
  my_theme +  geom_text(aes(label = ifelse(Count < 1, Count, sprintf("%.1f", Count))), family = "pt_sans_narrow",
                        position = position_dodge(width = 0.8), color = "black", size = 6, vjust = -0.3) +
  scale_y_continuous(limits = c(0, 4), breaks = seq(0, 4, 1), expand = c(0, 0)) +
  scale_fill_manual(values = c("Total"="#762a83","Novel"="#c2a5cf"))
p2

plot_grid <- p1 + p2 + plot_layout(axes = "collect", nrow = 1)
plot_grid

showtext_opts(dpi = 300)
ggsave(file=paste0("Fig1b.png"), width = 12, height = 3.5, bg="white", dpi = 300)
showtext_opts(dpi = 96)
showtext_auto()

