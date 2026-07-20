#ADMIXTURE PLOT 
#FIGURE : 2c 

setwd("/home/cuser48/Documents/GENOME_INDIA_PLOTS/GENOME_INDIA_CODES_TO_BE_SHARED_GIT/admixture_2c/")

#READING THE DATA :

data11_K5_modern <- read.csv("k5_admix_data.csv")

###################################################################################################################################
zz=as.data.frame(table(data11_K5_modern$pop_number))
colnames(zz)=c("pop","Freq")
zz=subset(zz,zz$Freq !=0)

zz2=as.data.frame(unique(data11_K5_modern$pop_number))
colnames(zz2)="pop"
zz2$numbers=rownames(zz2)
zz3=merge(zz2,zz,by="pop")
zz4=zz3[order(as.numeric(zz3$numbers)),]

z_K3=as.numeric(zz4$Freq)

z2_K3=vector()
sum=0

for (i in 1:length(z_K3)){
  sum=sum + z_K3[i]
  z2_K3[i]=sum 
}

z3_K3=unique(z2_K3)
z4_K3=c(0,z3_K3)

y_K3=vector()
for (i in 1:(length(z4_K3)-1)){
  y_K3[i]=z4_K3[i] + 0.5*(z4_K3[i+1]-z4_K3[i])
}

lab_num  = as.vector(unique(data11_K5_modern$pop_number))


###MAIN PLOT:

layout(matrix(c(1,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3),nrow=18))
par(mar=c(0,1,1,0))
par(xpd=F)

par(mar=c(0,2,1,0))

bar_positions <- barplot(rep(1,dim(data11_K5_modern)[1]), border = data11_K5_modern$Colour_LG_Tribe , yaxt="n" , xaxt = "n")
#axis(2,las=2,cex.axis=1.1,line=-2.5,font=2)
title(ylab = "LG", line = 0,cex.lab=1.8,font=2)
abline(v= bar_positions[z4_K3],col="white", cex =.1)

par(mar=c(0,2,1,0))

bar_positions <- barplot(rep(1,dim(data11_K5_modern)[1]), border = data11_K5_modern$Colour_Physiography , yaxt="n" , xaxt = "n")
#axis(2,las=2,cex.axis=1.1,line=-2.5,font=2)
title(ylab = "BG", line = 0,cex.lab=1.8,font=2)
abline(v=bar_positions[z4_K3],col="white", cex =.1)

par(mar=c(30,2,1,0))
# V1 = AAA (cluster_4), V2 = ATB (cluster_5), V3 = KURUMAN|PANIYAN (cluster_3) , V4 = ANI (cluster_1) , V5 = ASI (cluster_2)
original_colors <- c("steelblue", "palegreen3", "yellow3", "brown3", "navy")

# Function to create pastel colors
make_pastel <- function(color) {
  adjustcolor(color, alpha.f = 0.7) # Reduce opacity for pastel effect
}

# Generate pastel versions
pastel_colors <- sapply(original_colors, make_pastel)

# Plot to visualize the colors

barplot(t(as.matrix(data11_K5_modern[,1:5])), col=c( pastel_colors), xlab="", yaxt="n", border="NA",space=0.0, names.arg = rep("",dim(data11_K5_modern)[1]))
axis(1,at=z4_K3,labels=F)
Map(axis, side=1, line=0, at=y_K3, tick=F,labels=lab_num,las=2,srt =35,cex.axis=0.9,font=2)
axis(2,las=2,cex.axis=1.1,line=-2.5,font=2)
title(ylab = "Ancestry (K=5)", line = 0,cex.lab=1.8,font=2)
abline(v=z4_K3,col="white", cex =.1)
