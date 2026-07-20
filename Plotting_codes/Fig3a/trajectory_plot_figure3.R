set.seed(1234)
rm(list = ls())
library(latrend)
library(gridExtra)


NE_est_df <- read.csv("smcpp_Ne_estimates.csv")

NE_est_df <- NE_est_df[,c(1,2,3)]
NE_est_df[,c(2,3)] <- ceiling(NE_est_df[,c(2,3)])
colnames(NE_est_df) <- c("ID","Time","Ne")
NE_est_df <- NE_est_df[NE_est_df$Time>=85 & NE_est_df$Time <= 50000,]


Ne_df_gen <- data.frame(matrix(ncol = 3))
colnames(Ne_df_gen) <- c("ID","Time","Ne")

for (f in unique(NE_est_df$ID)){
  df = NE_est_df[NE_est_df$ID == f,]
  print(head(df))
  ids <- df[1,1]
  df <- df[,c(2,3)]
  Ne_df_fin <- data.frame(matrix(ncol = 2))
  colnames(Ne_df_fin) <- c("Time","Ne")
  est_pop_data <- df
  for (r1 in 1:(nrow(est_pop_data)-1)){
    #  print(r1)
    Ne_df <- est_pop_data[r1:(r1+1),]
    t <- Ne_df$Time
    Ne_t1 <- Ne_df$Ne[1]
    Ne_t2 <- Ne_df$Ne[2]
    alpha <- -log(Ne_t2/Ne_t1)/(diff(t))
    #  print(alpha)
    t <- seq(t[1]+1,t[2],by = 1)
    gen <- seq(1,length(t))
    if (alpha == 0){
      Ne_df_fin[t,] <- cbind(t,Ne_t1)
    }else{
      Ne_df_fin[t,] <- cbind(t,Ne_t1*exp(-alpha*gen))
    }
  }
  Ne_df_fin <- na.omit(Ne_df_fin)
  Ne_df_fin <- Ne_df_fin[Ne_df_fin$Time>=100 & Ne_df_fin$Time<=10000,]
  Ne_df_fin <- ceiling(Ne_df_fin)
  Ne_df_fin$ID <- ids
  Ne_df_fin <- Ne_df_fin[,c("ID","Time","Ne")]
  Ne_df_gen <- rbind(Ne_df_gen,Ne_df_fin)
}

Ne_df_gen <- na.omit(Ne_df_gen)
Ne_df_gen1 <- Ne_df_gen[Ne_df_gen$Time >= 300 & Ne_df_gen$Time<=2000,]
Ne_df_gen2 <- Ne_df_gen[Ne_df_gen$Time >= 300 & Ne_df_gen$Time<=10000,]

Ne_df_gen_temp <- Ne_df_gen[Ne_df_gen$Time >= 300 & Ne_df_gen$Time<=5000,]
Ne_df_gen_temp$Ne <- log10(Ne_df_gen_temp$Ne)

##############metadata df ###############333

metadata_df <- read.csv("metadata_Jan_2025.csv")
metadata_df <- metadata_df[,c("Code",
                              "LG_Tribe","Colour_LG_Tribe")]

metadata_df <- unique(metadata_df)

Ne_df_gen_temp <- merge(Ne_df_gen_temp,metadata_df, by.x = "ID",by.y = "Code")
Ne_df_gen_temp$LG_Tribe <- factor(Ne_df_gen_temp$LG_Tribe, levels = c("DR_NonTribe", "IE_NonTribe", "IE_Tribe",
                                              "TB_Tribe", "DR_Tribe", "AA_Tribe", "TB_NonTribe", "Unknown"))
colnames(Ne_df_gen_temp)
unique(Ne_df_gen_temp$ID)

head(Ne_df_gen_temp)

library(ggpubr)
library(ggplot2)
library(rlang)

plot_Ne <- function(df, group_col, color_col, na_label = "Unknown", na_color = "gray",outfile) {
  group_sym <- sym(group_col)
  color_sym <- sym(color_col)
  
  # Replace NA in group_col with 'Unknown' (or user-defined label)
  df <- df %>%
    mutate(
      !!group_sym := ifelse(is.na(!!group_sym), na_label, as.character(!!group_sym)),
      !!color_sym := ifelse(is.na(!!color_sym), na_color, as.character(!!color_sym))
    )
  
  # Build color mapping
  color_map <- unique(df[c(group_col, color_col)])
  color_vector <- setNames(color_map[[color_col]], color_map[[group_col]])
  color_vector <- color_vector[c("DR_NonTribe", "IE_NonTribe", "IE_Tribe", 
                                 "TB_Tribe", "DR_Tribe", "AA_Tribe", "TB_NonTribe", "Unknown")] 
  print(color_vector)
  p <- ggplot(df, aes(x = Time, y = Ne, group = ID)) +
    geom_line(aes(color = !!group_sym), size = 1, alpha = 0.8) +
    scale_color_manual(values = color_vector, name = group_col) +
    scale_x_continuous(breaks = seq(300, 5500, by = 500)) +
    theme_minimal() +
    labs(
      x = "Time (in generations)",
      y = "Effective population size (log10)",
    ) +
    theme(
      legend.position = "none",legend.key.size = 5,
      panel.grid.minor = element_blank())+
    theme_pubr(base_size =  20,base_family = "sans",legend = "none")

  ggsave(paste0(outfile,"_log10.png"),p,height = 6,width = 10, dpi = 300)
  return(p)
}
# Color by LG_Tribe using Colour_LG_Tribe
# appearance of overlapping line trajectories might differ from the provided figure due to
# use of Population code in plotting lines instead of population names

plot_Ne(
  df = Ne_df_gen_temp,
  group_col = "LG_Tribe",
  color_col = "Colour_LG_Tribe",
  outfile = "LG_Tribe_cluster_publication"
)
