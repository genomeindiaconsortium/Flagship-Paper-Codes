library(tidyverse)
library(vegan)

# Loading data ------------------------------------------------------------

#This input file is the supplementary table from the paper
metadata <- read.delim("./Supplementary_Table_S1.4.tsv", skip = 1) %>%
  select(Code, Numeric.Code, Linguistic.Group, Tribe) %>%
  mutate(Tribe = replace_values(Tribe,
                                "Yes" ~ "Tribe",
                                "No"  ~ "NonTribe",
                                 NA   ~ "CAO" )) %>%
  unite("group", Linguistic.Group, Tribe) %>%
  mutate(group = replace_values(group, "NA_CAO" ~ "CAO")) 


y <- read.delim("chr_Y_frequency.tsv") %>%
     column_to_rownames(var="X") 

mt <- read.delim("mtDNA_frequency.txt") %>% 
      pivot_wider(names_from = Haplogroup2, values_from = Proportion, values_fill = 0) %>% 
      column_to_rownames(var="Code")


# Calculating diversity ------------------------------------------------------


y_div <- as.data.frame(diversity(y, index = "shannon")) %>%
         rownames_to_column(var="Code") %>% 
         rename(y_diversity=2)


mt_div <- as.data.frame(diversity(mt, index = "shannon")) %>%
          rownames_to_column(var="Code") %>% 
          rename(mt_diversity=2)

df <- metadata %>% 
      full_join(mt_div, by ="Code") %>%
      full_join(y_div, by ="Code")


# Calculating Z-scores and saving ------------------------------------------------------

df <- df %>% 
  mutate(y_zscore = (y_diversity - mean(df$y_diversity)) / sd(df$y_diversity)) %>%
  mutate(mt_zscore = (mt_diversity - mean(df$mt_diversity)) / sd(df$mt_diversity)) %>%
  select(Code, Numeric.Code, group, mt_diversity, mt_zscore, y_diversity, y_zscore) %>%
  rename("Population" = 1,
         "Population Order" = 2,
         "Ethnolinguistic Group" = 3,
         "mtDNA Diversity" = 4,
         "mtDNA Diversity (Z-score)" = 5,
         "Y diversity" = 6,
         "Y Diversity (Z-score)" = 7)


write.table(df, "mtDNA_vs_chrY_AlphaBetaDiversity.tsv", quote = F, sep = "\t", row.names = F)

