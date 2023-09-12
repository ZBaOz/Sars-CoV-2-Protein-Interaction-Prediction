install.packages("xlsx")
library(xlsx)
setwd("/Veriler")

#READ POSITIVE PROTEINS FROM THE DATASET
pozitifProteinler <- read.xlsx(file = 'Proteinler.xlsx', 1, header=TRUE)
pozitifProteinler =as.character(as.vector(pozitifProteinler[,2]))

#READ ALL CANDIDATE PROTEINS FROM THE DATASET
adayProteinler <- read.xlsx(file = 'AdayProteinler.xlsx', 1, header=TRUE)
adayProteinler <- adayProteinler[,1]

#GET PROTEIN SEQUENCES FOR  POSITIVE PROTEINS AND CANDIDATE PROTEINS
install.packages("protr")
library("protr")
pozProtSeq<-getUniProt(pozitifProteinler)
adayProtSeq<-getUniProt(adayProteinler)

#LIBRARY DEFINITION FOR PAIRWISE ALIGNMENT
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("pairwiseAlignment")
library(Biostrings)
data(BLOSUM62)

install.packages("writexl")
library("writexl")

#identity values for positive and candidate proteins are keep in 'identity matrix'
identity = matrix(NA,nrow=length(pozitifProteinler),ncol=length(adayProteinler))
rownames(identity) = pozitifProteinler
colnames(identity) = adayProteinler
identity <- data.frame(identity)
dim(identity)

#The amino acid sequences of the proteins in the positive data set and the proteins in the candidate data set are drawn from the uniprot. 
#In the for loop, the sequence identity score of each positive protein with each candidate protein is calculated. The results are taken into the 'identity' matrix.
for (i in 332:length(pozitifProteinler))
{
  print(i)
  #s1 <- getUniProt(pozitifProteinler[i])
  s1=pozProtSeq[i]
  s1=unlist(s1,use.names=F)
  for (j in 1:length(adayProteinler))
  {
    #s2 <- getUniProt(adayProteinler[j])
    s2=adayProtSeq[j]
    s2=unlist(s2,use.names=F)
    palign1 <- pairwiseAlignment(s1, s2, substitutionMatrix = "BLOSUM62",type="global",gapOpening=8, gapExtension = 4)
    identity[i,j]=pid(palign1)
  }
  write_xlsx(identity,"/identity.xlsx")
}

identity2 <- read.xlsx(file = 'identity.xlsx', 1, header=TRUE)

#From positive interactions, data is read that captures how many human proteins interact with each Sars protein. 
#This is to generate as many negative interaction data as there are positive interactions with a sars protein.
library("readxl")
etkilesimSayi <- read.xlsx(file = 'Proteinler.xlsx', 7, header=TRUE)

#For example, there are 6 human proteins that interact with the SARS CoV 2 E protein. In this case, there will be 6 negative proteins. 
#The first 6 rows in the identity matrix contain the similarity of proteins interacting with cov2 e to candidate proteins. This 6 rows were taken into the protein matrix. 
#The aim is to determine the 6 candidate proteins with the lowest average similarity to cov2 e. The average of each column in the protein matrix is calculated. 
#Column names of min first 6 cells are labeled as non-interactive with sars cov2 e.
negatifProtein = matrix(NA,nrow=332,ncol=2)
indis=1
for (j in 1:length(etkilesimSayi))
{
  sayi=etkilesimSayi[1,j]
  protein=identity[1:sayi,]
  protein=protein[,2:5708]
  toplam=colSums(protein, na.rm = FALSE, dims = 1)
  ortalama=toplam/sayi
  # 
  # for(i in 1:(sayi*2))
  # {
  #   min=which.min(ortalama)
  #   ortalama[min]=10000
  # }
  
  for(i in 1:sayi)
  {
    min=which.min(ortalama)
    print(min)
    negatifProtein3[indis,1]=names(ortalama[min])
    negatifProtein3[indis,2]=ortalama[min]
    ortalama[min]=10000
    indis=indis+1
  }
  protein=NULL
}
#In the negative protein matrix, the first column is proteins, the 2nd column is average similarities.

negatifProtein <- data.frame(negatifProtein)
write_xlsx(negatifProtein,"negatifData.xlsx")



