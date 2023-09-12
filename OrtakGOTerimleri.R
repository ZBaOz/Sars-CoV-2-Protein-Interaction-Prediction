install.packages("xlsx")
library(xlsx)

install.packages("writexl")
library("writexl")

setwd("C:/Users/zeyne/Desktop/Kaggle Data/Proje2/Veriler")

#pozitif human proteinler için protein idleri al
pozHumProt <- read.xlsx(file = 'Proteinler.xlsx', 1, header=TRUE)
pozHumProt=as.character(as.vector(pozHumProt[,2]))

etkilesimSayi <- read.xlsx(file = 'Proteinler.xlsx', 7, header=TRUE)

covidProt <- read.xlsx(file = 'Proteinler.xlsx', 2, header=TRUE)
covidProt=as.character(as.vector(covidProt[1:27,3]))

predicted <- read.xlsx(file = 'C:/Users/zeyne/Desktop/Kaggle Data/Proje2/Makale1/Predicted Proteins.xlsx', 1, header=TRUE)



#uniprotr kütüphanesi protein ismine göre GO termleri çekmek için. 'GetProteinGOInfo' fonksiyonu ile yapýyor bunu
install.packages("UniprotR")
library("UniprotR")

#regex kullanmak için
#install.packages("stringr")
library("stringr")


#POZÝTÝF HUMAN PROTEÝNLER ÝLE ORTAK GO SAYILARI BULUNUR
ortakGOSayi=matrix(nrow=nrow(predicted),ncol=(length(etkilesimSayi)*3))

for (predIndis in 248:nrow(predicted))
#for (predIndis in 119:123)
{
  print(predIndis)
  k=1
  protIndis=1
  #bir covid proteinin etkileþimli olduðu tüm insan proteinlerinin 3 alanda GO terimleri çýkarýlýr.
  #dýþ döngü covid proteinlerini, iç döngü ise bir covid proteini ile etkileþimli insan proteinlerini dolaþýr.
  for (i in 1:length(etkilesimSayi))
  {
    sayi=etkilesimSayi[1,i]
    BPGOTerms=""
    CCGOTerms=""
    MFGOTerms=""
    #i. covid proteini ile etkileþimli insan proteinlerinin unique GO terimlerini çýkar
    for (j in 1:sayi)
    {
      #sýradaki proteinin GO terimlerini al
      GOObj<-GetProteinGOInfo(pozHumProt[protIndis])
      #BP tabanlý GO terimleri ekle
      x=str_extract_all(GOObj$Gene.ontology..biological.process,"GO:[0-9]+")
      x=unlist(x,use.names=F)
      BPGOTerms=c(BPGOTerms,x)
      
      #MF tabanlý GO terimleri ekle
      x2=str_extract_all(GOObj$Gene.ontology..molecular.function.,"GO:[0-9]+")
      x2=unlist(x2,use.names=F)
      MFGOTerms=c(MFGOTerms,x2)
      
      #CC tabanlý GO terimleri ekle
      x3=str_extract_all(GOObj$Gene.ontology..cellular.component.,"GO:[0-9]+")
      x3=unlist(x3,use.names=F)
      CCGOTerms=c(CCGOTerms,x3)
      
      protIndis=protIndis+1
    }
    BPGOTerms=unique(BPGOTerms)
    MFGOTerms=unique(MFGOTerms)
    CCGOTerms=unique(CCGOTerms)
    
    #Covid Proteinleri için GO terimler çýkarýlýr
    GOObjCov<-GetProteinGOInfo(covidProt[i])
    xBpCov=str_extract_all(GOObjCov$Gene.ontology..biological.process,"GO:[0-9]+")
    xBpCov=unlist(xBpCov,use.names=F)
    
    #MF tabanlý GO terimleri ekle
    xMfCov=str_extract_all(GOObjCov$Gene.ontology..molecular.function.,"GO:[0-9]+")
    xMfCov=unlist(xMfCov,use.names=F)
    
    #CC tabanlý GO terimleri ekle
    xCcCov=str_extract_all(GOObjCov$Gene.ontology..cellular.component.,"GO:[0-9]+")
    xCcCov=unlist(xCcCov,use.names=F)
    
    #pred protein için GO terimleri çýkarýlýr.
    #predIndis=6
    #GOObjH<-GetProteinGOInfo("A0A481SW80")
    GOObjH<-GetProteinGOInfo(predicted[predIndis,1])
    if (is.na(GOObjH$Gene.ontology..biological.process))
    {
      print("BP yok")
      ortakGOSayi[predIndis,k]=0
      covidGOSayi[predIndis,k]=0
    }
    else
    {
      #BP tabanlý GO terimlerini al, ortak BP terim sayýsýný matrise kaydet
      x=str_extract_all(GOObjH$Gene.ontology..biological.process,"GO:[0-9]+")
      x=unlist(x,use.names=F)
      a=intersect(BPGOTerms,x)
      ortakGOSayi[predIndis,k]=length(a)
      a=intersect(xBpCov,x)
      covidGOSayi[predIndis,k]=length(a)
    }
    if (is.na(GOObjH$Gene.ontology..molecular.function.))
    {
      print("MF yok")
      ortakGOSayi[predIndis,(k+1)]=0
      covidGOSayi[predIndis,(k+1)]=0
    }
    else
    {
      #MF tabanlý GO terimlerini al, ortak MF terim sayýsýný matrise kaydet
      x=str_extract_all(GOObjH$Gene.ontology..molecular.function.,"GO:[0-9]+")
      x=unlist(x,use.names=F)
      a=intersect(MFGOTerms,x)
      ortakGOSayi[predIndis,(k+1)]=length(a)
      a=intersect(xMfCov,x)
      covidGOSayi[predIndis,(k+1)]=length(a)
    }
    if (is.na(GOObjH$Gene.ontology..cellular.component.))
    {
      print("CC yok")
      ortakGOSayi[predIndis,(k+2)]=0
      covidGOSayi[predIndis,(k+2)]=0
    }
    else
    {
      #CC tabanlý GO terimlerini al, ortak CC terim sayýsýný matrise kaydet
      x=str_extract_all(GOObjH$Gene.ontology..cellular.component.,"GO:[0-9]+")
      x=unlist(x,use.names=F)
      a=intersect(CCGOTerms,x)
      ortakGOSayi[predIndis,(k+2)]=length(a)
      a=intersect(xCcCov,x)
      covidGOSayi[predIndis,(k+2)]=length(a)
    }
    k=k+3
  }
}

toplamGOSayi=matrix(nrow=nrow(predicted),ncol=length(etkilesimSayi))


#for (i in 1:nrow(predicted))
for (i in 1:285)
{
  indis=1
  for (k in 1:length(etkilesimSayi))
  {
    toplamGOSayi[i,k]=ortakGOSayi[i,indis]+ortakGOSayi[i,(indis+1)]+ortakGOSayi[i,(indis+2)]+covidGOSayi[i,indis]+covidGOSayi[i,(indis+1)]+covidGOSayi[i,(indis+2)]
    indis=indis+3
  }
}

toplamGOSayi2 <- data.frame(toplamGOSayi)
write_xlsx(toplamGOSayi2,"C:/Users/zeyne/Desktop/Kaggle Data/Proje2/Veriler/ToplamOrtakGOSayi2.xlsx")

