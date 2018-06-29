setwd("/Users/jleibowitz/Desktop/")
library('openxlsx')
library('gplots')
library('RColorBrewer')
source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/colAg.R")


data3<-read.xlsx('development.xlsx',rowNames=TRUE, sheet=1)



centered_data <- t(scale(t(data.matrix(theirdata)), scale = T, center = T))




test <- heatmap.2(centered_data, 
                  scale = "none", 
                  margins = c(5, 25), 
                  cexRow = 0.4, 
                  cexCol = 0.5, 
                  col = colorRampPalette(c("blue", "white", "red")),
                  xlab = "Samples", 
                  ylab = "Genes",
                  trace = "none",
                  keysize = 1,
                  key.title = NA,
                  Colv = FALSE,
                  Rowv = TRUE,
                  breaks = seq(from = -2, to = 2, by = 0.05),
                  main = "clusters",
                  key.xlab = "Down          Up",
                  density.info = "none", #remove histogram,
                  labRow = FALSE,
                  dendrogram='row'#, Colv=FALSE
)
cluster1<-labels(test$rowDendrogram[[2]][[2]])
write.xlsx(cluster1, file= "HCLcluster1.xlsx")
cluster2<-labels(test$rowDendrogram[[2]][[1]])
write.xlsx(cluster2, "HCLcluster2.xlsx")
cluster3<-labels(test$rowDendrogram[[1]][[2]])
write.xlsx(cluster3,"HCLcluster3.xlsx")
cluster4<-labels(test$rowDendrogram[[1]][[1]])
write.xlsx(cluster4,"HCLcluster4.xlsx")

centered_data2<-data.matrix(centered_data)
sorted <- centered_data2[match(rev(labels(test$rowDendrogram)), rownames(centered_data)), ]

# save the clusters and genes z-scores 
write.xlsx(sorted, 'CommongenesSODTREM2.xlsx',
           row.names = TRUE)




