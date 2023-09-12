install.packages("xlsx")
library(xlsx)
setwd("C:/Users/zeyne/Desktop/Covid-PPI Proje/Veriler")

#pozitif insan proteinleri dosyadan okunur
pozitifProteinler <- read.xlsx(file = 'Proteinler.xlsx', 1, header=TRUE)
pozitifProteinler =as.character(as.vector(pozitifProteinler[,2]))

#tüm aday proteinleri dosyadan okunur
adayProteinler <- read.xlsx(file = 'AdayProteinler.xlsx', 1, header=TRUE)
adayProteinler <- adayProteinler[,1]

#aday proteinler ile pozitif proteinlerin sekanslarý alýnýr.
install.packages("protr")
library("protr")
pozProtSeq<-getUniProt(pozitifProteinler)
adayProtSeq<-getUniProt(adayProteinler)

#pairwise alignment için kütüphane tanýmlarý
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("pairwiseAlignment")
library(Biostrings)
data(BLOSUM62)

install.packages("writexl")
library("writexl")

#pozitif veri kümesi ile aday proteinlerin identity deðerleri 'identity' matrisinde tutulur.
#matrisin satýr isimleri pozitif protein id'ler, sütun isimleri ise aday protein id ler olarak atandý
identity = matrix(NA,nrow=length(pozitifProteinler),ncol=length(adayProteinler))
rownames(identity) = pozitifProteinler
colnames(identity) = adayProteinler
identity <- data.frame(identity)
dim(identity)

#Pozitif veri kümesindeki proteinler ile aday veri kümesindeki proteinlerin amino asit sekanslarý
#uniprotdan çekilir. For döngüsünde, her bir pozitif proteinin her bir aday protein ile sequence identity
#skoru hesaplanýr. sonuçlar 'identity' matrisine alýnýr.
for (i in 332:length(pozitifProteinler))
{
  print(i)
  #s1 <- getUniProt(pozitifProteinler[183])
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
  write_xlsx(identity,"C:/Users/zeyne/Desktop/identity.xlsx")
}

identity2 <- read.xlsx(file = 'identity.xlsx', 1, header=TRUE)

#pozitif etkileþimlerden, her bir sars proteini ile kaç tane insan proteinin etkileþime girdiðini tutan veri
#okunur. bu bir sars proteini ile kaç tane pozitif etkileþim varsa o kadar negatif etkileþim verisi oluþturmak için
library("readxl")
etkilesimSayi <- read.xlsx(file = 'Proteinler.xlsx', 7, header=TRUE)

#örneðin sars cov 2 e ile etkileþimli 6 insan proteini var. Bu durumda 6 negatif protein bulunacak. identity
#matrisindeki ilk 6 satýr, cov2 e ile etikleþimli proteinlerin aday proteinler ile benzerliðini içerir.
#o 6 satýr protein matrisine alýndý. Amaç cov2 e ile en düþük ortalama benzerliðe sahip 6 aday proteini belirlemek
#protein matrsinde her bir sütunun ortalamasý hesaplanýr. min ilk 6 hücrenin sütun isimleri, sars cov2 e ile
#etkileþimsiz olarak etiketlenir.
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
#negatifProtein matrisinde ilk sütun proteinler 2. sütun ortalama benzerliklerdir.

negatifProtein <- data.frame(negatifProtein)
write_xlsx(negatifProtein,"C:/Users/zeyne/Desktop/negatifData.xlsx")



