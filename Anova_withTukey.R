my_data<-read.csv("/Users/jleibowitz/Desktop/APOEKO.csv",sep=",",na.strings=c("","NA"))
colnames(my_data)[1]<-"Groups"
Gene_names<-colnames(my_data)
Gene_names<-Gene_names[-1]

Anova_list <- list()
Anova_values <- data.frame()
Sig_genes <- data.frame()
post_hocs<-list()
Sig_pvalues <- data.frame()

for (i in 1:length(Gene_names)) {
  Anova_list[[i]] <- aov(my_data[,Gene_names[i]] ~ Groups, data=my_data)
  temp<-summary.aov(Anova_list[[i]])
  Anova_values[1,i] = temp[[1]][["Pr(>F)"]][[1]]
}

names(Anova_list)<-Gene_names

colnames(Anova_values)<-Gene_names
row.names(Anova_values)<-"pvalue"

t<-1
for (i in 1:length(Anova_values)){
  if(Anova_values[i] < 0.05){Sig_genes[1,t] <- names(Anova_values[i])
  Sig_pvalues[1,t] <- Anova_values[i]
  t=t+1
  }
}

for (i in 1:length(Sig_genes)){
  temp<-TukeyHSD(Anova_list[[Sig_genes[[i]]]])
  post_hocs[[i]]<-temp[[1]]
}

names(post_hocs)<-Sig_genes

numb_comps<-factorial(length(unique(my_data$Groups)))/(factorial(2)*factorial(length(unique(my_data$Groups))-2))
out<-data.frame(matrix(ncol = 1+numb_comps, nrow = length(Sig_genes)))
row.names(out)<-Sig_genes
temp<-combn(unique(my_data$Groups),2)
colnames(out) <- c("Anova p-value",t(rownames(post_hocs[[1]])))
out[,1]<-t(Sig_pvalues)

for(s in 1:length(post_hocs)){
  out[s,]<-c(out[s,1],t(post_hocs[[s]][,4]))
}

setwd("/Users/jleibowitz/Desktop")
write.csv(out,file=paste0("ANOVA_wTukey_",date(),".csv"))
