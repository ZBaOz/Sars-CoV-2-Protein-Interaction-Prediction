# -*- coding: utf-8 -*-
"""
Created on Fri Apr 15 20:25:47 2022

@author: zeyne
"""
import pandas as pd
import numpy as np

pip install xlrd
import xlrd
pip install openpyxl

#Word2Vec eğitim modeli oluşturmak için sekanslar yüklenir. 204.185 protein sekansı var. Uniprotdan indirdim
df = pd.read_csv(r'C:\Users\zeyne\Desktop\Covid-PPI Proje\Veriler\HomoSapSekans.csv')
sekanslar=df['V1']#sekans sütununun başlığı dosyada 'Entry'
#print (df['Entry'])
sekanslar.head
%reset
#POZİTİF VE NEGATİF PROTEİN SEKANSLARININ YÜKLENMESİ
tumData = pd.read_excel(r'C:\Users\zeyne\Desktop\Covid-PPI Proje\Veriler\VeriSeti3.xlsx', sheet_name='Sayfa1', engine='openpyxl')
PozProt=tumData['human sekans']
NegProt=tumData['negatif sekans']

#pozitif ve negatif proteinlerde Word2vec eğitim proteinlerine eklenir.204.849 protein oldu
tumData = pd.concat([sekanslar, PozProt], axis=0)
tumData = pd.concat([tumData, NegProt], axis=0)

tumData.head

#tokenization için kütüphane
import nltk 
WPT = nltk.WordPunctTokenizer()

from datetime import datetime
#çalışma süresi hesaplamak için get currenttime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

#Eğitim veri kümesi k-merlere bölünür
k=5#k aminoasitlik birimler birer kelime olacak
sekansSayi=len(tumData)
kmer_list = []#k-merler ile ifade edilmiş protein dokümanlarını tutacak bu liste
for i in range(sekansSayi):
    protein=tumData.iloc[i]
    start=0 # start ve end değişkenleri substringin başlangıç ve bitiş indisleri. 
    end=start+k 
    prot_k_mer=""#sıradaki proteinin k-merlerine bölünmüş halini bu değişken tutar. 
    k_merSayi=len(protein)-k+1
    j=0
    #her bir protein k harflik kelimelerden oluşan bir doküman olacak. k-mer ler birbirlerinden
    #birer boşlukla ayrılır tıpkı doğal dilde düz bir metin gibi.while döngüsünde her iterasyonda 
    #başlangıç ve bitiş 1 aminoasit kayar. yani k-meri alacak k birimlik pencere 1 aa kayıyordur
    while (j<k_merSayi):
        k_mer=protein[start:end:1]
        prot_k_mer=prot_k_mer+" "+k_mer
        start=start+1
        end=start+k
        j=j+1
    #'prot_k_mer' sıradaki proteinin k-merlerine ayrıştırılmış halini tutuyordur. başında fazladan bir
    #boşluk oluşuyor. 'lstrip' ile boşluk silinir.
    prot_k_mer=prot_k_mer.lstrip()
    #tokenize fonksiyonu ile tokenlerine ayrılır.
    tokens=WPT.tokenize(prot_k_mer)
    #tokenlerine ayrılmış protein k-mer listesine eklenir
    kmer_list.append(tokens)
len(kmer_list)    

from gensim.models import Word2Vec

#size=output vektörünün boyutu
#workers=paralelleştirebilmeye yarar. kaç core da paralel çalışacağı ('cython' yüklü değilse bir işe yaramaz)
#window=Maximum distance between the current and predicted word within a sentence.
#min_count= Ignores all words with total frequency lower than this.
#sg ({0, 1}– Training algorithm: 1 for skip-gram; otherwise CBOW. (doc2Vec'de distributed memory olduğu 
    #için burada cbow kullandım))
#hs ({0, 1}, optional) – If 1, hierarchical softmax will be used for model training. If 0, 
    #and negative is non-zero, negative sampling will be used.
#cbow_mean ({0, 1}, optional) – If 0, use the sum of the context word vectors. If 1, use the mean, 
    #only applies when cbow is used.
model5Mer = Word2Vec(kmer_list, min_count=1, vector_size=300, window=10, sg=0, negative=5, epochs=20,hs=0)
model3Mer.save(r"C:\Users\zeyne\Desktop\model3Mer.model")

word_vectors = model3Mer.wv
word_vectors.save(r"C:\Users\zeyne\Desktop\word2vec.wordvectors")
#print(model3Mer)
  
#POZİTİF PROTEİNLER
#oluşturulan vektörler ile proteinlerin temsil edilmesi için pozitif ve negatif proteinler de k-mer 
#lerine ayrılır.
k=5
sekansSayi=len(PozProt)
Pozkmer_list = []#k-merler ile ifade edilmiş protein dokümanlarını tutacak bu liste
for i in range(sekansSayi):
    protein=PozProt.iloc[i]
    start=0 # start ve end değişkenleri substringin başlangıç ve bitiş indisleri. 
    end=start+k 
    prot_k_mer=""#sıradaki proteinin k-merlerine bölünmüş halini bu değişken tutar. 
    k_merSayi=len(protein)-k+1
    j=0
    while (j<k_merSayi):
        k_mer=protein[start:end:1]
        prot_k_mer=prot_k_mer+" "+k_mer
        start=start+1
        end=start+k
        j=j+1
    prot_k_mer=prot_k_mer.lstrip()
    tokens=WPT.tokenize(prot_k_mer)
    Pozkmer_list.append(tokens)
    
#len(Pozkmer_list)   

data=pd.DataFrame()#poz ve neg veriyi tutacak olan değişken

#Dış döngü proteinlerde iç döngü ise bir proteinin tokenlerinde dolaşmayı sağlar. bir proteinindeki her
# tokenin vektörel karşılığı modelden alınır, 'sent_list' değişkenine eklenir. İç döngü her sonlandığında
#sent_list değişkeninde bir proteine ait tokenlerin vektörleri olacaktır. bu vektörlerin ortalaması 
#alınır. Ortalama vektörü proteinin Word2Vec vektörüdür, data değişkenine eklenir. 
for i in range (len(Pozkmer_list)):
    tokens=Pozkmer_list[i]
    sent_list=[]
    for word in tokens:
        wv2=model4Mer.wv[word]
        sent_list.append(wv2)
    vektor=np.mean(sent_list,axis=0)
    vektor=pd.DataFrame(vektor)
    vektor=np.transpose(vektor)
    data = data.append(vektor, ignore_index=True)

#NEGATİF PROTEİNLER
k=4
sekansSayi=len(NegProt)
Negkmer_list = []#k-merler ile ifade edilmiş protein dokümanlarını tutacak bu liste
for i in range(sekansSayi):
    protein=NegProt.iloc[i]
    start=0 # start ve end değişkenleri substringin başlangıç ve bitiş indisleri. 
    end=start+k 
    prot_k_mer=""#sıradaki proteinin k-merlerine bölünmüş halini bu değişken tutar. 
    k_merSayi=len(protein)-k+1
    j=0
    while (j<k_merSayi):
        k_mer=protein[start:end:1]
        prot_k_mer=prot_k_mer+" "+k_mer
        start=start+1
        end=start+k
        j=j+1
    prot_k_mer=prot_k_mer.lstrip()
    tokens=WPT.tokenize(prot_k_mer)
    Negkmer_list.append(tokens)


for i in range (len(Negkmer_list)):
    tokens=Negkmer_list[i]
    sent_list=[]
    for word in tokens:
        wv2=model4Mer.wv[word]
        sent_list.append(wv2)
    vektor=np.mean(sent_list,axis=0)
    vektor=pd.DataFrame(vektor)
    vektor=np.transpose(vektor)
    data = data.append(vektor, ignore_index=True)


#TEST VERİSİ
#test protein sekansları yükle
test = pd.read_csv(r'C:\Users\zeyne\Desktop\Covid-PPI Proje\Veriler\HomoSapSekans.csv')
test=test['V1']#sekans sütununun başlığı dosyada 'Entry'
print (test)
test.head
#test protein sekansları k-merlerine ayır
k=5
sekansSayi=len(test)
testkmer_list = []#k-merler ile ifade edilmiş protein dokümanlarını tutacak bu liste
for i in range(sekansSayi):
    protein=test.iloc[i]
    start=0 # start ve end değişkenleri substringin başlangıç ve bitiş indisleri. 
    end=start+k 
    prot_k_mer=""#sıradaki proteinin k-merlerine bölünmüş halini bu değişken tutar. 
    k_merSayi=len(protein)-k+1
    j=0
    while (j<k_merSayi):
        k_mer=protein[start:end:1]
        prot_k_mer=prot_k_mer+" "+k_mer
        start=start+1
        end=start+k
        j=j+1
    prot_k_mer=prot_k_mer.lstrip()
    tokens=WPT.tokenize(prot_k_mer)
    testkmer_list.append(tokens)

data=pd.DataFrame()

i=198277
#i=198277

#for i in range (len(testkmer_list)):
for i in range (198277,len(testkmer_list)):
    tokens=testkmer_list[i]
    sent_list=[]
    for word in tokens:
        wv2=model5Mer.wv[word]
        sent_list.append(wv2)
    vektor=np.mean(sent_list,axis=0)
    vektor=pd.DataFrame(vektor)
    vektor=np.transpose(vektor)
    data = data.append(vektor, ignore_index=True)

testkmer_list[194658]
data.iloc[194659]


data1=pd.DataFrame(data.iloc[0:100000,:])
data2=pd.DataFrame(data.iloc[100000:204186,:])
data3=pd.DataFrame(data.iloc[194658:204185,:])


#test verisi word2vec yazdırılması
writer = pd.ExcelWriter(r'C:\Users\zeyne\Desktop\negData3Word2Vec5Mer1.xlsx')
#write dataframe to excel
data1.to_excel(writer)
writer.save()

writer = pd.ExcelWriter(r'C:\Users\zeyne\Desktop\negData3Word2Vec5Mer2.xlsx')
#write dataframe to excel
data2.to_excel(writer)
writer.save()

writer = pd.ExcelWriter(r'C:\Users\zeyne\Desktop\negData3Word2Vec3Mer3.xlsx')
#write dataframe to excel
data3.to_excel(writer)
writer.save()

#compuation time için süre hesapla    
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)    

Label= pd.read_excel(r'C:\Users\zeyne\Desktop\Covid-PPI Proje\Veriler\Proteinler.xlsx', sheet_name='Sayfa3', header=0, engine='openpyxl')

data.insert(0,'etiket',Label)

writer = pd.ExcelWriter(r'C:\Users\zeyne\Desktop\word2Vec4Mer.xlsx')
# write dataframe to excel
data.to_excel(writer)
writer.save()



words = list(model3Mer.wv.vocab)
print(model4Mer.wv['MSSN'])



#Negatif örnekleme fikri, iyi bir modelin sahte sinyali gerçek sinyalden lojistik regresyon yoluyla 
#ayırt etmesi gerektiğini savunan gürültü karşılaştırmalı tahmin (benzer şekilde, üretici hasım ağları 
#gibi) kavramına dayanmaktadır. Ayrıca, negatif örnekleme hedefinin arkasındaki motivasyon, stokastik 
#gradyan inişine benzer: sahip olduğumuz binlerce gözlemin hepsini hesaba katarak her seferinde tüm 
#ağırlıkları değiştirmek yerine, sadece K tanesini kullanıyoruz ve hesaplamayı artırıyoruz. verimlilik 
#de önemli ölçüde (negatif numunelerin sayısına bağlıdır).


