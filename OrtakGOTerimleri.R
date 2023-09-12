install.packages("xlsx")
library(xlsx)

install.packages("writexl")
library("writexl")

setwd("/Veriler")

#Get protein ids for positive proteins.
pozHumProt <- read.xlsx(file = 'Proteinler.xlsx', 1, header=TRUE)
pozHumProt=as.character(as.vector(pozHumProt[,2]))

etkilesimSayi <- read.xlsx(file = 'Proteinler.xlsx', 7, header=TRUE)

covidProt <- read.xlsx(file = 'Proteinler.xlsx', 2, header=TRUE)
covidProt=as.character(as.vector(covidProt[1:27,3]))

predicted <- read.xlsx(file = 'Predicted Proteins.xlsx', 1, header=TRUE)



#'uniprotr' library is used finding common GO terms
install.packages("UniprotR")
library("UniprotR")

#library for regex
#install.packages("stringr")
library("stringr")


#FIND NUMBER OF COMMON GO TERMS WITH POSITIVE PROTEINS.
#THERE ARE 3 SUB ONTOLOGIES: BP, CC AND MF. THEREFORE NUMBER OF COLUMNS IS UP TO 3 TIMES THE NUMBER OF INTERACTIONS
ortakGOSayi=matrix(nrow=nrow(predicted),ncol=(length(etkilesimSayi)*3))

for (predIndis in 1:nrow(predicted))
{
  print(predIndis)
  k=1
  protIndis=1
  #GO terms are extracted for 3 sub-ontologies of all human proteins with which a covid protein interacts. The outer loop circulates covid proteins, 
  #and the inner loop circulates human proteins that interact with a covid protein.
  for (i in 1:length(etkilesimSayi))
  {
    sayi=etkilesimSayi[1,i]
    BPGOTerms=""
    CCGOTerms=""
    MFGOTerms=""
    #extract unique GO terms of human proteins interacting with i. covid protein
    for (j in 1:sayi)
    {
      #get GO terms of next proteins
      GOObj<-GetProteinGOInfo(pozHumProt[protIndis])
      #BP GO terms
      x=str_extract_all(GOObj$Gene.ontology..biological.process,"GO:[0-9]+")
      x=unlist(x,use.names=F)
      BPGOTerms=c(BPGOTerms,x)
      
      #MF GO terms
      x2=str_extract_all(GOObj$Gene.ontology..molecular.function.,"GO:[0-9]+")
      x2=unlist(x2,use.names=F)
      MFGOTerms=c(MFGOTerms,x2)
      
      #CC GO terms
      x3=str_extract_all(GOObj$Gene.ontology..cellular.component.,"GO:[0-9]+")
      x3=unlist(x3,use.names=F)
      CCGOTerms=c(CCGOTerms,x3)
      
      protIndis=protIndis+1
    }
    BPGOTerms=unique(BPGOTerms)
    MFGOTerms=unique(MFGOTerms)
    CCGOTerms=unique(CCGOTerms)
    
    #Extract GO terms fo Covid Proteins
    GOObjCov<-GetProteinGOInfo(covidProt[i])
    #BP GO terms
    xBpCov=str_extract_all(GOObjCov$Gene.ontology..biological.process,"GO:[0-9]+")
    xBpCov=unlist(xBpCov,use.names=F)
    
    #MF GO terms
    xMfCov=str_extract_all(GOObjCov$Gene.ontology..molecular.function.,"GO:[0-9]+")
    xMfCov=unlist(xMfCov,use.names=F)
    
    #CC GO terms
    xCcCov=str_extract_all(GOObjCov$Gene.ontology..cellular.component.,"GO:[0-9]+")
    xCcCov=unlist(xCcCov,use.names=F)
    
    #Extract GO terms for predicted proteinsd.
    #predIndis=6
   
    GOObjH<-GetProteinGOInfo(predicted[predIndis,1])
    if (is.na(GOObjH$Gene.ontology..biological.process))
    {
      print("No BP term")
      ortakGOSayi[predIndis,k]=0
      covidGOSayi[predIndis,k]=0
    }
    else
    {
      #Get BP GO terms. Save number of common BP terms to the matrix.
      x=str_extract_all(GOObjH$Gene.ontology..biological.process,"GO:[0-9]+")
      x=unlist(x,use.names=F)
      a=intersect(BPGOTerms,x)
      ortakGOSayi[predIndis,k]=length(a)
      a=intersect(xBpCov,x)
      covidGOSayi[predIndis,k]=length(a)
    }
    if (is.na(GOObjH$Gene.ontology..molecular.function.))
    {
      print("no MF")
      ortakGOSayi[predIndis,(k+1)]=0
      covidGOSayi[predIndis,(k+1)]=0
    }
    else
    {
      #Get MF GO terms. Save number of common MF terms to the matrix.
      x=str_extract_all(GOObjH$Gene.ontology..molecular.function.,"GO:[0-9]+")
      x=unlist(x,use.names=F)
      a=intersect(MFGOTerms,x)
      ortakGOSayi[predIndis,(k+1)]=length(a)
      a=intersect(xMfCov,x)
      covidGOSayi[predIndis,(k+1)]=length(a)
    }
    if (is.na(GOObjH$Gene.ontology..cellular.component.))
    {
      print("no CC")
      ortakGOSayi[predIndis,(k+2)]=0
      covidGOSayi[predIndis,(k+2)]=0
    }
    else
    {
      #Get CC GO terms. Save number of common CC terms to the matrix.
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


for (i in 1:nrow(predicted))
{
  indis=1
  for (k in 1:length(etkilesimSayi))
  {
    toplamGOSayi[i,k]=ortakGOSayi[i,indis]+ortakGOSayi[i,(indis+1)]+ortakGOSayi[i,(indis+2)]+covidGOSayi[i,indis]+covidGOSayi[i,(indis+1)]+covidGOSayi[i,(indis+2)]
    indis=indis+3
  }
}

toplamGOSayi2 <- data.frame(toplamGOSayi)
write_xlsx(toplamGOSayi2,"ToplamOrtakGOSayi2.xlsx")

