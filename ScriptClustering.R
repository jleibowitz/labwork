setwd("/Users/jleibowitz/Desktop/Monocytes__>microglia")
library('openxlsx')
library('gplots')
library('RColorBrewer')
source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/colAg.R")

# get the data and fold change
data <- read.xlsx('FDR_APOE34remove_outlier.xlsx', rowNames = TRUE)

# I calculated the log2 fold change manually and saved it to this file 
FC <- read.xlsx('FDR_APOE34remove_outlier.xlsx', sheet = 2)

countDFrpkm_mean <- colAg(myMA=data, group=c(rep(1, 2), rep(2, 3), rep(3, 3), rep(4, 3)), myfct=mean)
countDFrpkm_mean<-data.frame(countDFrpkm_mean)
write.xlsx(countDFrpkm_mean, "test.xlsx", row.names=TRUE) 

# get the differential genes - FC > 2 
APOE3_P3.P14 <- FC[FC[,2] >= 0.5 | FC[,2] <= -0.5,] 
APOE4_P3.P14 <- FC[FC[,3] >= 0.5 | FC[,3] <= -0.5,] 
APOE3_P3.APOE4_P3 <- FC[FC[,4] >= 0.5 | FC[,4] <= -0.5,] 
APOE3_P14.APOE4_P14 <- FC[FC[,5] >= 0.5 | FC[,5] <= -0.5,] 
APOE3_P3.APOE4_P14 <- FC[FC[,6] >= 0.5 | FC[,6] <= -0.5,] 
APOE3_P14.APOE4_P3 <- FC[FC[,7] >= 0.5 | FC[,7] <= -0.5,] 

# get all the unique genes 

sigdiff <- unique(c(APOE3_P3.P14$Name, APOE4_P3.P14$Name, APOE3_P3.APOE4_P3$Name, APOE3_P14.APOE4_P14$Name, APOE3_P3.APOE4_P14$Name, APOE3_P14.APOE4_P3$Name))

class(sigdiff)

# 750 diff genes
final <- data[sigdiff, ] # each replicate 
#final <- countDFrpkm_mean[sigdiff, ] # get the averages 
final<-read.xlsx("mgmo.xlsx",sheet=14,rowNames=T)

# calculate the z-score and generate heatmap 
centered_data <- t(scale(t(data.matrix(final)), scale = T, center = T))
test <- heatmap.2(centered_data, 
                  scale = "none", 
                  #margins = c(7, 10), 
                  cexRow = 0.4, 
                  cexCol = 0.5, 
                  col = colorRampPalette(rev(brewer.pal(n = 8, name = "RdBu"))),
                  xlab = "Samples", 
                  ylab = "Genes",
                  trace = "none",
                  keysize = 1,
                  key.title = NA,
                  Colv = TRUE,
                  Rowv = TRUE,
                  breaks = seq(from = -2, to = 2, by = 0.15),
                  main = "clusters",
                  key.xlab = "Down          Up",
                  density.info = "none", #remove histogram,
                  labRow = FALSE,
                  dendrogram='both'#, Colv=FALSE
)


# perform k-means
set.seed(12345)
centers <- 20
final <- data.frame(centered_data)

# Assign genes to a cluster and sort clusters by 1-10
km <- kmeans(final, centers = centers, iter.max=10)
m.kmeans<- cbind(final, cluster = km$cluster)
order <- m.kmeans[order(m.kmeans$cluster),]
colnames(order) <- c('APOE3 P3', 'APOE4 P3', 'APOE3 P14', 'APOE4 P14', 'cluster')

# arrange clusters based on size (highest --> low )
len <- c()
for(i in 1:centers){
  len <- c(len, length(which(order$cluster == i)))
}

len2 <- order(-len)

sort.final <- data.frame()
for ( x in len2){
  sort.final <- rbind(sort.final, order[which(order$cluster == x),])
}

# assign colors to clusters -- can change 
colours <-  colorRampPalette(brewer.pal(11,"Spectral"))(centers)[sort.final$cluster]

# generate heatmap with clusters 
test <- heatmap.2(data.matrix(sort.final[, -c(ncol(sort.final))]), 
                  scale = "none", 
                  margins = c(10, 30), 
                  cexRow = 0.4, 
                  cexCol = 0.5, 
                  col = colorRampPalette(rev(brewer.pal(n = 8, name = "RdBu"))),
                  xlab = "Samples", 
                  ylab = "Genes",
                  trace = "none",
                  keysize = 1,
                  key.title = NA,
                  Colv = FALSE,
                  Rowv = FALSE,
                  #Rowv = TRUE,
                  breaks = seq(from = -2, to = 2, by = 0.15),
                  main = "clusters",
                  key.xlab = "Down          Up",
                  density.info = "none", #remove histogram,
                  labRow = FALSE,
                  RowSideColors=colours,
                  dendrogram='none'#, Colv=FALSE
)

legend("bottomright",      
       legend = unique(sort.final$cluster),
       col = unique(colours), 
       lty= 1,             
       lwd = 5,           
       cex=.7
)


# save the clusters and genes z-scores 
write.xlsx(sort.final, 'mycluster.xlsx',
           row.names = TRUE)




