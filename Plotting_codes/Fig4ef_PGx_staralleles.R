



pharmacogenomics_input <- as.data.frame(read.csv("Pharmacogenomics_GI9768_14genes.phenotype.mod.tsv", header = TRUE, sep = "\t"))
head(pharmacogenomics_input)

pharmacogenomics_input %>%
  filter(Phenotype != "ultrarapid_metabolizer") %>%
  ggplot(aes(x = Gene, fill = Phenotype)) +
  geom_bar(position = "fill", stat = "count", width = 0.9) +
  theme_classic(20) +
  theme(axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 20, color = "black", face = "italic"),
        plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(limits = rev) +
  scale_fill_manual(name = "Metabolizer", labels = c("Intermediate", "Normal", "Poor", "Rapid", "Unknown"),
                    values = c("#FF8200", "#06629E", "#A51C28", "#0C7C59", "#666666")) +
  labs(y = "Proportion", title = "") +
  coord_flip()

showtext_opts(dpi = 300)
ggsave("./4e_staralleles_all.png", width = 6.5, height = 6, dpi = 300, bg="white")
showtext_opts(dpi = 96)
showtext_auto()


group_order <- c("AA_T", "DR_NT", "DR_T", "IE_NT", "IE_T", "TB_NT", "TB_T", "CAO")

selectedGenes_input <- subset(
  pharmacogenomics_input,
  Gene %in% c("VKORC1", "UGT1A1", "SLCO1B1", "CYP3A5", "CYP2C19", "CYP2B6"))
selectedGenes_input$Phenotype <- factor(selectedGenes_input$Phenotype,
                                        levels = c("intermediate_metabolizer", "normal_metabolizer",
                                                   "poor_metabolizer", "rapid_metabolizer",
                                                   "ultrarapid_metabolizer", "unknown_metabolizer"))
selectedGenes_input$Group <- factor(selectedGenes_input$Group,
                                    levels = rev(levels(selectedGenes_input$Group)))
#plot
selectedGenes_plot1 <- ggplot(selectedGenes_input, aes(x=fct_relevel(Group, rev(group_order)), fill=Phenotype)) +
  geom_bar(position = "fill", stat = "count", width = 0.8) +
  facet_wrap(~Gene, nrow=2) +
  scale_fill_manual(name = "Metabolizer", labels = c("Intermediate", "Normal", "Poor", "Rapid",
                               "Unknown"),
                    values = c("#FF8200", "#06629E", "#A51C28", "#0C7C59", "#666666")) +
  theme_classic(20) +
  guides(fill = guide_legend(ncol=1)) +
  labs(y = "Proportion") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=18, family="Arial", color = "black"),
        strip.text.x = element_text(size=18, family="Arial", color = "black", face = "italic"),
        axis.title.y = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_rect(fill = "lightgrey", color = "black", linewidth = 1),
        legend.position = "right",
        legend.title = element_text(size=20, family="Arial"),
        legend.text = element_text(size=18, family="Arial"),
        legend.box.spacing = unit(0.5, "cm"),
        legend.key.size = unit(0.8, "cm")) +
  coord_flip()

selectedGenes_plot1

showtext_opts(dpi = 300)
ggsave("./4e_staralleles_subset.png", width = 11, height = 8, dpi = 300, bg="white")
showtext_opts(dpi = 96)
showtext_auto()

selectedGenes_plot2 <- selectedGenes_input %>%
  filter(Phenotype != "ultrarapid_metabolizer") %>%
  ggplot(aes(x=fct_relevel(Group, rev(group_order)), fill=Phenotype)) +
  geom_bar(position = "fill", stat = "count", width = 0.8) +
  facet_wrap(~Gene, nrow=1) +
  scale_fill_manual(name = "Metabolizer", labels = c("Intermediate", "Normal", "Poor", "Rapid",
                                                     "Unknown"),
                    values = c("#FF8200", "#06629E", "#A51C28", "#0C7C59", "#666666")) +
  theme_classic(20) +
  # guides(fill = guide_legend(nrow=1)) +
  labs(y = "Proportion") +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=18, family="Arial", color = "black"),
        strip.text.x = element_text(size=18, family="Arial", color = "black", face = "italic"),
        axis.title.y = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_rect(fill = "lightgrey", color = "black", linewidth = 1),
        legend.position = "bottom",
        legend.title = element_text(size=20, family="Arial"),
        legend.text = element_text(size=18, family="Arial"),
        legend.box.spacing = unit(0.5, "cm"),
        legend.key.size = unit(1, "cm")) +
  coord_flip()

selectedGenes_plot2

showtext_opts(dpi = 300)
ggsave("./4e_staralleles_subset.png", width = 15, height = 6, dpi = 300, bg="white")
showtext_opts(dpi = 96)
showtext_auto()
