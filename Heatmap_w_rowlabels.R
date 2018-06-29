setwd("/Users/jleibowitz/Desktop")
library('openxlsx')
library('gplots')
library('RColorBrewer')
source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/colAg.R")


data<-read.xlsx('APOEKOM.xlsx',rowNames=TRUE)


centered_data <- t(scale(t(data.matrix(data)), scale = T, center = T))

incl_gene_names<-vector(mode="character")
incl_gene_names<-c("Tmem119", "Olfml3","Siglech","Fcrls","P2ry12","Arg1","Arg2","Il1b","Il1rn") #INPUT GENE NAMES TO INCLUDE ON HEAT MAP
rowlabels<-vector(mode="character", length=nrow(data))

for (i in 1:length(incl_gene_names)){
  rowlabels[grep(incl_gene_names[i],rownames(data))]<-incl_gene_names[i]
  
}


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
                  labRow = rowlabels,
                  dendrogram='row'#, Colv=FALSE
)


centered_data2<-data.matrix(centered_data)
sorted <- centered_data2[match(rev(labels(test$rowDendrogram)), rownames(centered_data)), ]

# save the clusters and genes z-scores 
write.xlsx(sorted, 'CommongenesSODTREM2.xlsx',
           row.names = TRUE)




