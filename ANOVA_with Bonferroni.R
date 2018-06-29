my_data<-read.csv("/Users/jleibowitz/Desktop/Test.csv",sep=",",na.strings=c("","NA"))
colnames(my_data)[1]<-"Groups"
Gene_names<-colnames(my_data)
Gene_names<-Gene_names[-1]
i <- 1
Anova_values <- data.frame()
Sig_genes <- data.frame()
Sig_pvalues <- data.frame()

while (i < length(Gene_names)+1) {
  Anova_temp <- aov(my_data[,Gene_names[i]] ~ Groups, data=my_data)
  temp<-summary.aov(Anova_temp)
  Anova_values[1,i] = temp[[1]][["Pr(>F)"]][[1]]
  i = i+1
  
}

colnames(Anova_values)<-Gene_names
row.names(Anova_values)<-"pvalue"
i<-1
t<-1
while (i < length(Anova_values)+1){
  if(Anova_values[i] < 0.05){Sig_genes[1,t] <- names(Anova_values[i])
  Sig_pvalues[1,t] <- Anova_values[i]
  t=t+1
  }
  i = i+1
}

post_hocs<-list()

for(i in 1:length(Sig_genes)){
  temp<-pairwise.t.test(my_data[,Sig_genes[1,i]], my_data$Groups, p.adj="bonferroni")  
  post_hocs[[i]]<-temp[[3]]
}
names(post_hocs)<-Sig_genes

numb_comps<-factorial(length(unique(my_data$Groups)))/(factorial(2)*factorial(length(unique(my_data$Groups))-2))
out<-data.frame(matrix(ncol = 1+numb_comps, nrow = length(Sig_genes)))
row.names(out)<-Sig_genes
temp<-combn(unique(my_data$Groups),2)

temp2<-"Anova p-value"
for(i in 1:ncol(post_hocs[[1]])){
  for(t in 1:nrow(post_hocs[[1]])){
    if(is.na(post_hocs[[1]][t,i])){next}
    temp2<-c(temp2, paste0(colnames(post_hocs[[1]])[i],"/",row.names(post_hocs[[1]])[t]))
  }
}


colnames(out) <- temp2
out[,1]<-t(Sig_pvalues)

for(s in 1:length(post_hocs)){
  p<-2
  for(i in 1:ncol(post_hocs[[s]])){
    for(t in 1:nrow(post_hocs[[s]])){
      if(is.na(post_hocs[[s]][t,i])){next}
      out[s,p]<-post_hocs[[s]][t,i]
      p=p+1
    }
  }
}

setwd("/Users/jleibowitz/Desktop")
write.csv(out,file=paste0("ANOVA_wPosthoc_",date(),".csv"))
