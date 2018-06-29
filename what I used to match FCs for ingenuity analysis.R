library('openxlsx')
setwd("/Users/jleibowitz/Desktop")
data<-read.xlsx('Apoe diff regulated genes.xlsx',rowNames=FALSE,colNames = TRUE)

for (i in 1:length(data$Unique.Up.in.males)){
  x<-grep(data$Unique.Up.in.males[i],data$`Males-Up`)
  data$`FC-Unique.Up.Males`[i]<-data$`FC-Males.Up`[x]
  
}

for (i in 1:length(data$Unique.Up.in.males)){
  x<-grep(data$Unique.up.in.females[i],data$`Females-up`)
  data$`FC-Unique.Up.in.Females`[i]<-data$`FC-Females.Up`[x]
  
}

for (i in 1:length(data$Unique.Up.in.males)){
  x<-grep(data$Unique.Down.in.males[i],data$`Males-Down`)
  data$`FC-Unique.Down.in.males`[i]<-data$`FC-Males.Down`[x]
  
}

for (i in 1:length(data$Unique.Up.in.males)){
  x<-grep(data$Unique.down.in.females[i],data$`Females-Down`)
  data$`FC-Unique.down.in.Females`[i]<-data$`FC-Females.Down`[x]
  
}

for (i in 1:length(data$Unique.Up.in.males)){
  x<-grep(data$Common.Down[i],data$`Males-Down`)
  data$`FC-Common.Down`[i]<-data$`FC-Males.Down`[x]
  
}

for (i in 1:length(data$Unique.Up.in.males)){
  x<-grep(data$Common.Down[i],data$`Females-Down`)
  data$`FC-Common.Down`[i]<-data$`FC-Females.Down`[x]
  
}

for (i in 1:length(data$Unique.Up.in.males)){
  x<-grep(data$Common.Up[i],data$`Males-Up`)
  data$`FC-Common.Up`[i]<-data$`FC-Males.Up`[x]
  
}

for (i in 1:length(data$Unique.Up.in.males)){
  x<-grep(data$Common.Up[i],data$`Females-up`)
  data$`FC-Common.Up`[i]<-data$`FC-Females.Up`[x]
  
}

write.table(data$`FC-Common.Up`, "clipboard.csv",sep=",")

write.xlsx(data, 'APOE_diff_regulated_genes',
           row.names = TRUE)
