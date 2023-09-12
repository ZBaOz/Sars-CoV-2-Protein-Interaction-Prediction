install.packages("protr")
library("protr")

install.packages("xlsx")
library(xlsx)

setwd("Veriler/")

data <- read.xlsx(file = 'VeriSeti.xlsx', 1, header=TRUE)
#select pozitive proteins from dataset
pozProt =as.character(as.vector(data[,8]))
#select negative proteins from dataset
negProt =as.character(as.vector(data[,9]))

#EXTRACTING AAC FEATURES. 

#AAC features for positive proteins
sekans=unlist(pozProt[1],use.names=F)
AACOzellik=extractAAC(sekans)
for (i in c(1:332))
{
  sekans=unlist(pozProt[i],use.names=F)
  y=extractAAC(sekans)
  AACOzellik = rbind(AACOzellik, y)
}

#Extracting AAC features for Negative Proteins
for (i in c(1:332))
{
  sekans=unlist(negProt[i],use.names=F)
  y=extractAAC(sekans)
  AACOzellik = rbind(AACOzellik, y)
}
write.xlsx(AACOzellik2, '/AACozellik2.xlsx', sheetName = "Sheet1", col.names = TRUE, row.names = FALSE, append = FALSE)

data <- read.xlsx(file = 'VeriSeti2.xlsx', 1, header=TRUE)
negProt =as.character(as.vector(data[,6]))

sekans=unlist(negProt[1],use.names=F)
AACOzellik2=extractAAC(sekans)
for (i in c(2:332))
{
  sekans=unlist(negProt[i],use.names=F)
  y=extractAAC(sekans)
  AACOzellik2 = rbind(AACOzellik2, y)
}

data <- read.xlsx(file = 'VeriSeti3.xlsx', 1, header=TRUE)
negProt =as.character(as.vector(data[,6]))

sekans=unlist(negProt[1],use.names=F)
AACOzellik2=extractAAC(sekans)
for (i in c(2:332))
{
  sekans=unlist(negProt[i],use.names=F)
  y=extractAAC(sekans)
  AACOzellik2 = rbind(AACOzellik2, y)
}


#EXTRACTING PAAC FEATURES. 'extractPAAC' FUNCTION HAS SOME PARAMETERS."Hydrophobicity", "Hydrophilicity" AND "SideChainMass" ARE DEFAULT PARAMETERS. 
#ALL PARAMETERS WERE USED WITH THEIR DEFAULT VALUES. "Lamda" PARAMETER DETERMINE THE FEATURE SIZE. BY SETTING THIS VALUE TO 30 SO 50 FEATURES DATASET WAS COMPOSED. 

data <- read.xlsx(file = 'VeriSeti.xlsx', 1, header=TRUE)
pozProt =as.character(as.vector(data[,8]))

negProt =as.character(as.vector(data[,9]))

#PAAC features for positive proteins
sekans=unlist(pozProt[1],use.names=F)
PozPAACOzellik=extractPAAC(sekans, props = c("Hydrophobicity", "Hydrophilicity", "SideChainMass"), lambda = 30, w = 0.05, customprops = NULL)
for (i in c(2:332))
{
  sekans=unlist(pozProt[i],use.names=F)
  y=extractPAAC(sekans, props = c("Hydrophobicity", "Hydrophilicity","SideChainMass"), lambda = 30, w = 0.05, customprops = NULL)
  PozPAACOzellik = rbind(PozPAACOzellik, y)
}

#PAAC features for negative proteins. since there are three negative datasets, same procedure repeated 3 times.
neg1PAAC=PozPAACOzellik

for (i in c(1:332))
{
  sekans=unlist(negProt[i],use.names=F)
  y=extractPAAC(sekans, props = c("Hydrophobicity", "Hydrophilicity","SideChainMass"), lambda = 30, w = 0.05, customprops = NULL)
  neg1PAAC = rbind(neg1PAAC, y)
}

write.xlsx(neg1PAAC, 'PAACozellik.xlsx', sheetName = "Sheet1", col.names = TRUE, row.names = FALSE, append = FALSE)

data <- read.xlsx(file = 'VeriSeti2.xlsx', 1, header=TRUE)
negProt =as.character(as.vector(data[,6]))

neg2PAAC=PozPAACOzellik

for (i in c(1:332))
{
  sekans=unlist(negProt[i],use.names=F)
  y=extractPAAC(sekans, props = c("Hydrophobicity", "Hydrophilicity","SideChainMass"), lambda = 30, w = 0.05, customprops = NULL)
  neg2PAAC = rbind(neg2PAAC, y)
}

write.xlsx(neg2PAAC, 'PAACozellik2.xlsx', sheetName = "Sheet1", col.names = TRUE, row.names = FALSE, append = FALSE)

data <- read.xlsx(file = 'VeriSeti3.xlsx', 1, header=TRUE)
negProt =as.character(as.vector(data[,6]))

neg3PAAC=PozPAACOzellik

for (i in c(1:332))
{
  sekans=unlist(negProt[i],use.names=F)
  y=extractPAAC(sekans, props = c("Hydrophobicity", "Hydrophilicity", "SideChainMass"), lambda = 30, w = 0.05, customprops = NULL)
  neg3PAAC = rbind(neg3PAAC, y)
}

write.xlsx(neg3PAAC, 'PAACozellik2.xlsx', sheetName = "Sheet1", col.names = TRUE, row.names = FALSE, append = FALSE)



