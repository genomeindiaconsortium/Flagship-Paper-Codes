rm(list = ls())
admx_df <- read.csv("admix_file_k5.txt",sep = " ")
admx_df <- admx_df[,2:7]

colnames(admx_df)
# match the columns while running
# [1] "ID"  "V1"  "V2"  "V3"  "V4"  "V5" 
library(vegan)  # for diversity() function

# Select only the columns with ancestry proportions
ancestry_matrix <- admx_df[, c("V1", "V2", "V3", "V4", "V5")]

# Calculate Shannon diversity for each row (individual)
admx_df$shannon_alpha <- diversity(ancestry_matrix, index = "shannon")
colnames(admx_df)[1] <- "IID" 

##########33
metadata_df <- read.csv("GI_metadata.csv")
metadata_df1 <- metadata_df[,c(2,5,6,19,21,30)]
colnames(metadata_df1)

#[1] "IID"               "Code"              "Population_Order"          "Admixture_Cluster"
#[6] "Tribe_Shape"       "LG_Tribe"  

admx_df1 <- merge(admx_df,metadata_df1,by = "IID")

head(admx_df1[1:5,])

# Calculate per-population mean alpha and beta
alpha_summary <- aggregate(shannon_alpha ~ Code, data = admx_df1, FUN = median)

# For beta: calculate pairwise distances and average them per population
library(vegan)
ancestry_matrix <- admx_df1[, c("V1", "V2", "V3", "V4", "V5")]
dist_matrix <- as.matrix(vegdist(ancestry_matrix, method = "bray"))
admx_df1$Code <- as.factor(admx_df1$Code)

# Compute average pairwise distance (beta) per population
beta_summary <- do.call(rbind, lapply(levels(admx_df1$Code), function(pop) {
  ids <- which(admx_df1$Code == pop)
  if (length(ids) > 1) {
    beta_val <- median(dist_matrix[ids, ids][lower.tri(dist_matrix[ids, ids])])
  } else {
    beta_val <- NA  # Not enough individuals
  }
  data.frame(Code = pop, beta = beta_val)
}))

# Combine alpha and beta summaries
div_summary <- merge(alpha_summary, beta_summary, by = "Code")

# Standardize
div_summary$alpha_z <- scale(div_summary$shannon_alpha)[,1]
div_summary$beta_z  <- scale(div_summary$beta)[,1]


div_summary$class <- with(div_summary, ifelse(alpha_z < 0 & beta_z < 0, "Low α,Low β",
                                              ifelse(alpha_z < 0 & beta_z > 0, "Low α,High β",
                                                     ifelse(alpha_z > 0 & beta_z < 0, "High α,Low β",
                                                            "High α,High β"))))


colnames(unique(metadata_df[,c(5,6,19,21,25)]))

#> colnames(unique(metadata_df[,c(5,6,19,21,25)]))
#[1] "Code"                     "Population_Order"         "Admixture_Cluster"       
#[4] "Tribe_Shape"              "Colour_Admixture_Cluster"
 
div_summary <- merge(div_summary, unique(metadata_df[,c(5,6,19,21,25)]), by = "Code")
write.table(div_summary,"alpha_beta_admixture_diversity_paper.txt",col.names = T,row.names = F,quote = F,sep = "\t")
