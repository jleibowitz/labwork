setwd("/Users/jleibowitz/Desktop")
library('openxlsx')
library('gplots')
library('RColorBrewer')
library('pca3d')
library('rgl')
my_data<-read.xlsx("heatmap.xlsx",sheet = 2)
data<-my_data[,-1]
pca<-prcomp(data, scale.=TRUE)

plotPCA <- function(x, nGroup) {
  n <- ncol(x) 
  if(!(n %in% c(2,3))) { # check if 2d or 3d
    stop("x must have either 2 or 3 columns")
  }
  
  fit <- hclust(dist(x), method="complete") # cluster
  groups <- cutree(fit, k=nGroup)
  
  if(n == 3) { # 3d plot
    plot3d(x, col=groups, type="s", size=1, axes=F)
    axes3d(edges=c("x--", "y--", "z"), lwd=3, axes.len=2, labels=FALSE)
    grid3d("x")
    grid3d("y")
    grid3d("z")
  } else { # 2d plot
    maxes <- apply(abs(x), 2, max)
    rangeX <- c(-maxes[1], maxes[1])
    rangeY <- c(-maxes[2], maxes[2])
    plot(x, col=groups, pch=19, xlab=colnames(x)[1], ylab=colnames(x)[2], xlim=rangeX, ylim=rangeY)
    lines(c(0,0), rangeX*2)
    lines(rangeY*2, c(0,0))
  }
}


fit<-hclust(dist(pca$x[,1:3]), method="complete")
groups<-cutree(fit, k=4)

plotPCA(pca$x[,1:3],5)

#plotting using pca3d
gr<-factor(my_data[,1])
pca3d(pca, group=gr, legend= "topleft", fancy=TRUE)
pca2d(pca, group=gr, legend="topleft")
#####################################