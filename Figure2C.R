setwd("/Users/jleibowitz/Desktop")
library('openxlsx')
library('gplots')
library('RColorBrewer')
source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/colAg.R")

# get the data and fold change
data <- read.xlsx('APOEKO.xlsx', rowNames = TRUE)


# calculate the z-score and generate heatmap 
centered_data <- t(scale(t(data.matrix(data)), scale = T, center = T))


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
                  Rowv = FALSE,
                  breaks = seq(from = -2, to = 2, by = 0.05),
                  main = "clusters",
                  key.xlab = "Down          Up",
                  density.info = "none", #remove histogram,
                  labRow = FALSE,
                  dendrogram='none'#, Colv=FALSE
)


centered_data2<-data.matrix(centered_data)
sorted <- centered_data2[match(rev(labels(test$rowDendrogram)), rownames(centered_data)), ]

# save the clusters and genes z-scores 
write.xlsx(sorted, 'CommongenesSODTREM2.xlsx',
           row.names = TRUE)




