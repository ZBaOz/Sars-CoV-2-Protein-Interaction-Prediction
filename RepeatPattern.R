install.packages("xlsx")
library(xlsx)

#library for regex
install.packages("stringr")
library("stringr")

install.packages("writexl")
library("writexl")

install.packages("protr")
library("protr")

#aminoacid vector
aa <- c('A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V')

setwd("Veriler/")


#REPEATED PATTERNS IN EACH PROTEINS
#SINCE THERE ARE 20 AMINO ACIDS, EACH PROTEIN IS REPRESENTED BY A 20-DIMENSIONAL VECTOR.

#Read positive proteins from the dataset
pozProt <- read.xlsx(file = 'Proteinler.xlsx', 3, header=TRUE)
pozProt=as.character(as.vector(pozProt[,2]))

#Repeat Pattern features for positive protein
RepPatPoz1=matrix(nrow=length(pozProt),ncol=20)
colnames(RepPatPoz1) = aa
RepPatPoz1 <- data.frame(RepPatPoz1)

for (i in 1:length(pozProt))
{
  protein=pozProt[i]
  for (j in 1:20)
  {
    tmp <- c(aa[j],"+")
    str=paste(tmp, collapse="")
    yer=str_extract_all(protein,str)
    yer=unlist(yer,use.names=F)
    toplam=0
    if (length(yer)>0)
    {
      for (k in 1:length(yer))
      {
        sayi=nchar(yer[k])
        toplam=toplam+(sayi*sayi)
      }
    }
    RepPatPoz1[i,j]=toplam
  }
}
write_xlsx(RepPatPoz1,"RepPattern.xlsx")

#read covid proteins from the dataset
covidProt <- read.xlsx(file = 'Proteinler.xlsx', 2, header=TRUE)
covidProt=as.character(as.vector(covidProt[,2]))

#Repeated aminoacids for each covid protein
RepPatCovid=matrix(nrow=length(covidProt),ncol=20)
colnames(RepPatCovid) = aa
RepPatCovid <- data.frame(RepPatCovid)

for (i in 1:length(covidProt))
{
  protein=covidProt[i]
  for (j in 1:20)
  {
    tmp <- c(aa[j],"+")
    str=paste(tmp, collapse="")
    yer=str_extract_all(protein,str)
    yer=unlist(yer,use.names=F)
    toplam=0
    if (length(yer)>0)
    {
      for (k in 1:length(yer))
      {
        sayi=nchar(yer[k])
        toplam=toplam+(sayi*sayi)
      }
    }
    RepPatCovid[i,j]=toplam
  }
}
write_xlsx(RepPatCovid,"RepPatternCov.xlsx")


#Read Negative dataset from the dataset.
#Since there are 3 negative datasets this process should repeat for each negative dataset
negProt <- read.xlsx(file = 'negatifData.xlsx', 1, header=TRUE)
negProt=as.character(as.vector(negProt[,1]))

#Repeated aminoacids for each negative protein
RepPatNeg=matrix(nrow=length(negProt),ncol=20)
colnames(RepPatNeg) = aa
RepPatNeg <- data.frame(RepPatNeg)

for (i in 1:length(negProt))
{
  proteinID=negProt[i]
  protein<-getUniProt(proteinID)
  for (j in 1:20)
  {
    tmp <- c(aa[j],"+")
    str=paste(tmp, collapse="")
    yer=str_extract_all(protein,str)
    yer=unlist(yer,use.names=F)
    toplam=0
    if (length(yer)>0)
    {
      for (k in 1:length(yer))
      {
        sayi=nchar(yer[k])
        toplam=toplam+(sayi*sayi)
      }
    }
    RepPatNeg[i,j]=toplam
  }
}
write_xlsx(RepPatNeg,"RepPatternNeg.xlsx")

Sys.time()



#EACH PROTEIN SPLITTED T 5 EQUAL SIZE SUB SEQUENCE. REPEAT PATTERNS ARE EXTRACTED FOR EACH SUB SEQUENCE 
#THERE ARE 0 AMINO ACIDS AND 5 SUB SEQUENCE FOR EACH PROTEIN. THEREFORE EACH PROTEIN IS REPRESENTED WITH 5X20 SIZE VECTOR

sutun=c(aa,aa,aa,aa,aa)


#POSITIVE PROTEINS
pozProt <- read.xlsx(file = 'Proteinler.xlsx', 3, header=TRUE)
pozProt=as.character(as.vector(pozProt[,2]))

RepPatPozParcali=matrix(nrow=length(pozProt),ncol=100)
colnames(RepPatPozParcali) = sutun
RepPatPozParcali <- data.frame(RepPatPozParcali)

for (i in 1:length(pozProt))
{
  protein=pozProt[i]
  uz=nchar(protein)
  parcaSize=ceiling(uz/5)
  bas=1;
  for (m in 1:4)
  {
    son=parcaSize*m
    parca[m]=substr(protein, bas, son)
    bas=son+1
  }
  parca[5]=substr(protein, bas, uz)
  
  stIndis=1
  for (m in 1:5)
  {
    for (j in 1:20)
    {
      tmp <- c(aa[j],"+")
      str=paste(tmp, collapse="")
      yer=str_extract_all(parca[m],str)
      yer=unlist(yer,use.names=F)
      toplam=0
      if (length(yer)>0)
      {
        for (k in 1:length(yer))
        {
          sayi=nchar(yer[k])
          toplam=toplam+(sayi*sayi)
        }
      }
      RepPatPozParcali[i,stIndis]=toplam
      stIndis=stIndis+1
    }
  }
}

write_xlsx(RepPatPozParcali,"RepPattern.xlsx")

#COVID PROTEINS
covidProt <- read.xlsx(file = 'Proteinler.xlsx', 2, header=TRUE)
covidProt=as.character(as.vector(covidProt[,2]))


RepPatCovidParcali=matrix(nrow=length(covidProt),ncol=100)
colnames(RepPatCovidParcali) = sutun
RepPatCovidParcali <- data.frame(RepPatCovidParcali)

for (i in 1:length(covidProt))
{
  protein=covidProt[i]
  uz=nchar(protein)
  parcaSize=ceiling(uz/5)
  bas=1;
  for (m in 1:4)
  {
    son=parcaSize*m
    parca[m]=substr(protein, bas, son)
    bas=son+1
  }
  parca[5]=substr(protein, bas, uz)
  
  stIndis=1
  for (m in 1:5)
  {
    for (j in 1:20)
    {
      tmp <- c(aa[j],"+")
      str=paste(tmp, collapse="")
      yer=str_extract_all(parca[m],str)
      yer=unlist(yer,use.names=F)
      toplam=0
      if (length(yer)>0)
      {
        for (k in 1:length(yer))
        {
          sayi=nchar(yer[k])
          toplam=toplam+(sayi*sayi)
        }
      }
      RepPatCovidParcali[i,stIndis]=toplam
      stIndis=stIndis+1
    }
  }
}

write_xlsx(RepPatCovidParcali,"RepPatternCov.xlsx")


#NEGATIVE PROTEINS
negProt <- read.xlsx(file = 'negatifData.xlsx', 1, header=TRUE)
negProt=as.character(as.vector(negProt[,1]))

RepPatNegParcali=matrix(nrow=length(negProt),ncol=100)
colnames(RepPatNegParcali) = sutun
RepPatNegParcali <- data.frame(RepPatNegParcali)

for (i in 1:length(negProt))
{
  proteinID=negProt[i]
  protein<-getUniProt(proteinID)
  uz=nchar(protein)
  parcaSize=ceiling(uz/5)
  bas=1;
  for (m in 1:4)
  {
    son=parcaSize*m
    parca[m]=substr(protein, bas, son)
    bas=son+1
  }
  parca[5]=substr(protein, bas, uz)
  
  stIndis=1
  for (m in 1:5)
  {
    for (j in 1:20)
    {
      tmp <- c(aa[j],"+")
      str=paste(tmp, collapse="")
      yer=str_extract_all(parca[m],str)
      yer=unlist(yer,use.names=F)
      toplam=0
      if (length(yer)>0)
      {
        for (k in 1:length(yer))
        {
          sayi=nchar(yer[k])
          toplam=toplam+(sayi*sayi)
        }
      }
      RepPatNegParcali[i,stIndis]=toplam
      stIndis=stIndis+1
    }
  }
}
Sys.time()

write_xlsx(RepPatNegParcali,"RepPatternNeg.xlsx")

