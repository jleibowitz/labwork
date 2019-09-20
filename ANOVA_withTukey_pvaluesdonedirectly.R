setwd("/users/jleibowitz/Dropbox (Partners Healthcare)/Oleg and Jeff/6 - Monocytes-->microglia in EAE/Heatmapanalyses")
data<-read.xlsx("mgmo.xlsx",sheet=2,colNames = FALSE,rowNames = FALSE)
data[,2]<-NULL
rownames(data)<-data[,1]
data[,1]<-NULL
rownames(data)[1]<-"Groups"

my_data<-t(data)

Gene_names<-colnames(my_data)
Gene_names<-Gene_names[-1]
Gene_names<-as.vector(Gene_names)

my_data<-data.frame(my_data)
Anova_list <- list()
Anova_values <- data.frame()
Sig_genes <- data.frame()
post_hocs<-list()
Sig_pvalues <- data.frame()

for (i in 1:length(Gene_names)){
  as.numeric(as.character(my_data$Gene_names[i]))
}

for (i in 1:length(Gene_names)) {
  Anova_list[[i]] <- aov(as.numeric(as.character(my_data[,Gene_names[i]])) ~ Groups, data=my_data)
  temp<-summary.aov(Anova_list[[i]])
  Anova_values[1,i] = temp[[1]][["Pr(>F)"]][[1]]
}

names(Anova_list)<-Gene_names

colnames(Anova_values)<-Gene_names
row.names(Anova_values)<-"pvalue"

final<-data.frame(t(Anova_values))
write.xlsx(final,"pvalues.xlsx")







