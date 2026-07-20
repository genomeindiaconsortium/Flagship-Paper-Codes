set.seed(1234)
library(latrend)
library(gridExtra)

colors_ls = c("#000000", "#3cb44b" ,"#ffe119", "#4363d8", "#f58231", "#911eb4", 
              "#ffff00", "#f032e6", "#bcf60c", "#b03060", "#3f000f", "#a57c00","#fabebe", 
              "#008080", "#ff4500", "#9a6324", "#b8860b", "#800000", "#aaffc3", "#808000", 
              "#ffd8b1", "#000075", "#808080", "#dc143c","#6b8e23","#00ffff","#008080",
              "#c71585","#0000cd","#deb887","#a020f0","#adff2f","#2f4f4f","#e9967a",
              "#8b4513","#191970","#8b0000","#778899","#bc8f8f","#228b22","#4682b4",
              "#cd5c5c","#8b008b","#1e90ff","#90be6d","#dc143c")

#loading smc++ data for all population groups
NE_est_df <- read.csv("smcpp_Ne_estimates.csv")
NE_est_df <- NE_est_df[,c(1,2,3)]
NE_est_df[,c(2,3)] <- ceiling(NE_est_df[,c(2,3)])
colnames(NE_est_df) <- c("ID","Time","Ne")
NE_est_df <- NE_est_df[NE_est_df$Time>=85 & NE_est_df$Time <= 50000,]

Ne_df_gen <- data.frame(matrix(ncol = 3))
colnames(Ne_df_gen) <- c("ID","Time","Ne")

# function to interpolate data within points (t, t+1) based on change in 
# population size from timepoint t to t+1
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

# subset data to do clustering on data upto 2000 and 10000 generations

Ne_df_gen1 <- Ne_df_gen[Ne_df_gen$Time >= 300 & Ne_df_gen$Time<=2000,]
Ne_df_gen2 <- Ne_df_gen[Ne_df_gen$Time >= 300 & Ne_df_gen$Time<=10000,]

#function to do clustering on longitudinal Ne trajectories

gckmbase <- lcMethodGCKM(
  formula = Ne ~ (Time|ID),
  id = "ID",
  time = "Time",control = lmerControl(optimizer ="Nelder_Mead"),iter.max = 10000)

# function to select best K
choose_K <- function(Ne_df,method_obj,outfile){
  gckm_methods <- lcMethods(method = method_obj, nClusters = 1:20)
  
  models <- latrendBatch(gckm_methods, data = Ne_df)
  metric_df <- metric(models, c("WMAE", "BIC"))
  #plotMetric(models, c("WMAE", "BIC"))
  p1 <- plotMetric(models, c("WMAE"))+
    labs(x = "K")
  ggsave(paste0(outfile,"_WMAE.png"),plot = p1,height = 6,width = 6,dpi = 300)
  write.table(metric_df,paste0(outfile,"_metric.txt"),col.names = T,row.names = F,quote = F,sep ="\t")
}

set.seed(1234)
choose_K(Ne_df=Ne_df_gen1,method_obj=gckmbase,outfile="Ne_300_2000_K")
choose_K(Ne_df=Ne_df_gen2,method_obj=gckmbase,outfile="Ne_300_10000_K")

#function to perform clustering at a specific choice of K and save the files
fitted_trajectory_K <- function(K,Ne_df,method_obj,outfile){
  K_colors <- sample(colors_ls,K)
  gckm <- latrend(method_obj, data = Ne_df,nClusters = K)
  #plot(gckm)
  p1 <- plotClusterTrajectories(gckm,linewidth = 1) + 
    scale_color_manual(values = K_colors)+
    labs(x = "Time(in generations)")
  p1
  ggsave(paste0(outfile,"_fittedTrajectory.png"),plot = p1,height = 6,width = 8,dpi = 300)
  df_fin <- fittedTrajectories(gckm)
  df_fin$ID <- as.factor(df_fin$ID)
  colnames(df_fin) <- c("ID","Time","Ne_fitted","Cluster")
  Ne_df <- merge(Ne_df,df_fin,by = c("ID","Time"))
  write.table(Ne_df,paste0(outfile,"_fittedTrajectory.txt"),col.names = T,row.names = F,quote = F,sep ="\t")
}

set.seed(1234)
for (k in 4:8) {
  fitted_trajectory_K(K = k,Ne_df = Ne_df_gen1,method_obj = gckmbase,outfile = paste0("test_Ne_300_2000_K_",k))
  fitted_trajectory_K(K = k,Ne_df = Ne_df_gen2,method_obj = gckmbase,outfile = paste0("Ne_300_10000_K_",k))
}
