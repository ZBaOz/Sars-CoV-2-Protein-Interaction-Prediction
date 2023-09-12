install.packages("protr")
library("protr")

install.packages("xlsx")
library(xlsx)
setwd("/Veriler/")


#POSITIVE, NEGATIVE AND COVID PROTEIN SEQUENCE ARE READ FROM DATASET TO SEPERATE VECTORS.
data <- read.xlsx(file = 'VeriSeti3.xlsx', 1, header=TRUE)
pozProt =as.character(as.vector(data[,2]))

negProt =as.character(as.vector(data[,3]))

covidProt=as.character(as.vector(data[,1]))


Sys.time()

#CONJOINT TRIAD FOR COVID PROTEINS
#'extractCTriad' FUNCTION FROM 'protr' LIBRARY WAS USED. AFTER FEATURE EXTRACTION PROCESS A 343-FEATURE DATASET IS OBTAINED.  
sekans=unlist(covidProt[1],use.names=F)
CTriadOzellik=extractCTriad(sekans)
for (i in c(2:332))
{
  sekans=unlist(covidProt[i],use.names=F)
  y=extractCTriad(sekans)
  CTriadOzellik = rbind(CTriadOzellik, y)
}
write.xlsx(CTriadOzellik, 'CTriadOzellikPoz.xlsx', sheetName = "Sheet1", col.names = TRUE, row.names = FALSE, append = FALSE)


#CONJOINT TRIAD FOR NEGATIVE PROTEINS
sekans=unlist(negProt[1],use.names=F)
CTriadOzellik=extractCTriad(sekans)
for (i in c(2:332))
{
  sekans=unlist(negProt[i],use.names=F)
  y=extractCTriad(sekans)
  CTriadOzellik = rbind(CTriadOzellik, y)
}


#CONJOINT TRIAD FOR POSITIVE PROTEINS
sekans=unlist(pozProt[1],use.names=F)
CTriadOzellik=extractCTriad(sekans)
for (i in c(2:332))
{
  sekans=unlist(pozProt[i],use.names=F)
  y=extractCTriad(sekans)
  CTriadOzellik = rbind(CTriadOzellik, y)
}

Sys.time()

