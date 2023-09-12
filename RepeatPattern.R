install.packages("xlsx")
library(xlsx)



#regex kullanýmý için
install.packages("stringr")
library("stringr")

install.packages("writexl")
library("writexl")

#negatif protinlerin sekanslarýný çekmek için
install.packages("protr")
library("protr")

#dogadaki 20 aminoasit
aa <- c('A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V')

setwd("C:/Users/zeyne/Desktop/Covid-PPI Proje/Veriler/")


#HER BÝR PROTEÝNDEKÝ TEKRAR EDEN AA PATTERN

#Pozitif Proteinler
#pozitif protein sekanslarý dosyada 3. sheette 2. sütunda
pozProt <- read.xlsx(file = 'Proteinler.xlsx', 3, header=TRUE)
pozProt=as.character(as.vector(pozProt[,2]))

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
write_xlsx(RepPatPoz1,"C:/Users/zeyne/Desktop/RepPattern.xlsx")

#Covid Proteinleri
covidProt <- read.xlsx(file = 'Proteinler.xlsx', 2, header=TRUE)
covidProt=as.character(as.vector(covidProt[,2]))

#her bir proteindeki tekrar eden single aminoasitler
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
write_xlsx(RepPatCovid,"C:/Users/zeyne/Desktop/RepPatternCov.xlsx")


#Negatif Proteinler
negProt <- read.xlsx(file = 'negatifData3.xlsx', 1, header=TRUE)
negProt=as.character(as.vector(negProt[,1]))

#her bir proteindeki tekrar eden single aminoasitler
RepPatNeg=matrix(nrow=length(negProt),ncol=20)
colnames(RepPatNeg) = aa
RepPatNeg <- data.frame(RepPatNeg)

for (i in 284:length(negProt))
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
write_xlsx(RepPatNeg,"C:/Users/zeyne/Desktop/RepPatternNeg.xlsx")

Sys.time()

#parca=matrix()

#HER BÝR PROTEÝN 5 PARÇAYA BÖLÜNDÜ HER BÝR PARÇADAKÝ TEKRAR EDEN AA PATTERN

sutun=c(aa,aa,aa,aa,aa)


#Pozitif Proteinler
#pozitif protein sekanslarý dosyada 3. sheette 2. sütunda
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

write_xlsx(RepPatPozParcali,"C:/Users/zeyne/Desktop/RepPattern.xlsx")

#Covid Proteinleri
covidProt <- read.xlsx(file = 'Proteinler.xlsx', 2, header=TRUE)
covidProt=as.character(as.vector(covidProt[,2]))

covidProt=covidProt[-c(32,33)]
covidProt[31]

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

write_xlsx(RepPatCovidParcali,"C:/Users/zeyne/Desktop/RepPatternCov.xlsx")


#Negatif Proteinler
negProt <- read.xlsx(file = 'negatifData3.xlsx', 1, header=TRUE)
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

write_xlsx(RepPatNegParcali,"C:/Users/zeyne/Desktop/RepPatternNeg.xlsx")

covidProt <- covidProt[-c(32,33)]
