# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 19:39:17 2022

@author: zeyne
"""

import pandas as pd
import numpy as np
#Doc2Vec eğitim modeli oluşturmak için sekanslar yüklenir. 204.185 protein sekansı var. Uniprotdan indirdim
df = pd.read_csv(r'C:\Users\zeyne\Desktop\Covid-PPI Proje\Veriler\HomoSapProteins.csv')
sekanslar=df['Entry']#sekans sütununun başlığı dosyada 'entry'
print (df['Entry'])
df.head
%reset
#POZİTİF PROTEİNLERİN YÜKLENMESİ
dfPozitif = pd.read_excel(r'C:\Users\zeyne\Desktop\Covid-PPI Proje\Veriler\Proteinler.xlsx', sheet_name='Human Proteins')
PozProt=dfPozitif['Protein Sequence']

#NEGATİF PROTEİNLERİN YÜKLENMESİ
dfNegatif = pd.read_excel(r'C:\Users\zeyne\Desktop\Covid-PPI Proje\Veriler\negatifData3.xlsx', sheet_name='Sayfa1')
NegProt=dfNegatif['sekans']

#pozitif ve negatif proteinlerde doc2vec eğitim proteinlerine eklenir.204.849 protein oldu
sekanslar = pd.concat([sekanslar, PozProt], axis=0)
sekanslar = pd.concat([sekanslar, NegProt], axis=0)

#train test bölmek için rastgele indis numarası üretir
import random
indis=np.arange(start=0, stop=332, step=1)
PozRowNumTest=random.sample(range(0, 332), 67)
PozRowNumTrain=np.delete(indis, PozRowNumTest)
NegRowNumTest=random.sample(range(0, 332), 67)
NegRowNumTrain=np.delete(indis, NegRowNumTest)


from datetime import datetime
#çalışma süresi hesaplamak için get currenttime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

#DOC2VEC EĞİTMEK İÇİN KULLANILACAK PROTEİNLERDEN K-MER ÜRETME
k=3#k aminoasitlik birimler birer kelime olacak
sekansSayi=len(sekanslar)
kmer_list = []#k-merler ile ifade edilmiş protein dokümanlarını tutacak bu liste
for i in range(sekansSayi):
    protein=sekanslar.iloc[i]
    start=0 # start ve end değişkenleri substringin başlangıç ve bitiş indisleri. 
    end=start+k
    prot_k_mer="[" #sıradaki proteinin k-merlerine bölünmüş halini bu değişken tutar. 
    k_merSayi=len(protein)-k+1
    j=0
    #her bir protein k harflik kelimelerden oluşan bir doküman olacak. Ancak doc2vec fonksiyonuna bunu verebilmek
    #için '["ILR", "IGN"]' şeklinde [] parantezleri arasında, her biri virgülle ayrılmış ve çift tırnak arasında
    #bulunacak kmer ler şeklinde düzenlemeyi sağlar while döngüsü. Her iterasyonda başlangıç ve bitiş 1 aminoasit
    # kayar
    while (j<k_merSayi-1):
        k_mer=protein[start:end:1]
        prot_k_mer=prot_k_mer+"'"+k_mer+"',"+" "
        start=start+1
        end=start+k
        j=j+1
    #proteinin son k-mer i fazladan bir virgül oluşmasın diye döngü dışında eklenir
    k_mer=protein[start:end:1]
    prot_k_mer=prot_k_mer+"'"+k_mer+"'"
    prot_k_mer=prot_k_mer+"]"
    #ilgili protein için k-mer dokümanı listeye eklemir
    kmer_list.append(prot_k_mer)

len(kmer_list)
sekanslar.iloc[204183]
 
#DOC2VEC EĞİTİP MODEL OLUŞTURMA
from gensim.test.utils import common_texts
from gensim.models.doc2vec import Doc2Vec, TaggedDocument

del documents

documents = [TaggedDocument(doc, [i]) for i, doc in enumerate(kmer_list)]
len(documents)
#dm eğitim algoritma seçimi. 1 ise distributed memory dir.
#vector_size= çıktı vektörü kaç boyutlu olacak
#window=The maximum distance between the current and predicted word within a sentence.
#Ignores all words with total frequency lower than this.
modelDM3mer = Doc2Vec(documents, vector_size=300, window=10, min_count=1, workers=4, dm=1, epochs=20, negative=5)


now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

#modelin kaydedilmesi
modelDM3mer.save("C:/Users/zeyne/Desktop/modelDM7mer.mod")
vector = model4mer.infer_vector(["ILR", "IGN"])
model = Doc2Vec.load("./doc2vecmodel.mod")


#EĞİTİLEN DOC2VEC MODELİ İLE NEGATİF PROTEİNLERİN VEKTÖRLERE DÖNÜŞTÜRÜLMESİ
#Pozitif proteinlerin doküman vektörleri elde edilecek. öncelikle proteinlerin k-merlerine ayrılması gerek. proteinler
# k ardışık aa dan oluşan kelimeler ile, aralarında birer boşlukla ifade edilir
k=5#k aminoasitlik birimler birer kelime olacak
sekansSayi=len(PozProt)
kmer_listPoz = []#k-merler ile ifade edilmiş pozitif protein dokümanlarını tutacak bu liste
vektorPoz=[]
vektorPoz=pd.DataFrame(vektorPoz)
for i in range(sekansSayi):
    protein=PozProt[i]
    start=0 # start ve end değişkenleri substringin başlangıç ve bitiş indisleri. 
    end=start+k
    prot_k_mer="" #sıradaki proteinin k-merlerine bölünmüş halini bu değişken tutar. 
    k_merSayi=len(protein)-k+1
    j=0
    while (j<k_merSayi):
        k_mer=protein[start:end:1]
        prot_k_mer=prot_k_mer+" "+k_mer
        start=start+1
        end=start+k
        j=j+1
    li = list(prot_k_mer.split(" "))
    vector = modelDM5mer.infer_vector(li)
    vector=pd.DataFrame(vector)
    vector=np.transpose(vector)
    vektorPoz = vektorPoz.append(vector, ignore_index=True)
    #ilgili protein için k-mer dokümanı listeye eklemir
    kmer_listPoz.append(prot_k_mer)


#PozTest=vektorPoz.iloc[PozRowNumTest,:] 
#PozTrain=vektorPoz.iloc[PozRowNumTrain,:] 
#PozTest.shape

vektorPoz=pd.DataFrame(vektorPoz)
writer = pd.ExcelWriter(r'C:\Users\zeyne\Desktop\Pozitif.xlsx')
#write dataframe to excel
vektorPoz.to_excel(writer)
writer.save()
    

#EĞİTİLEN DOC2VEC MODELİ İLE NEGATİF PROTEİNLERİN VEKTÖRLERE DÖNÜŞTÜRÜLMESİ
k=5#k aminoasitlik birimler birer kelime olacak
sekansSayi=len(NegProt)
kmer_listNeg = []#k-merler ile ifade edilmiş pozitif protein dokümanlarını tutacak bu liste
vektorNeg=[]
vektorNeg=pd.DataFrame(vektorNeg)
for i in range(sekansSayi):
    protein=NegProt[i]
    start=0 # start ve end değişkenleri substringin başlangıç ve bitiş indisleri. 
    end=start+k
    prot_k_mer="" #sıradaki proteinin k-merlerine bölünmüş halini bu değişken tutar. 
    k_merSayi=len(protein)-k+1
    j=0
    while (j<k_merSayi):
        k_mer=protein[start:end:1]
        prot_k_mer=prot_k_mer+" "+k_mer
        start=start+1
        end=start+k
        j=j+1
    li = list(prot_k_mer.split(" "))
    vector = modelDM4mer.infer_vector(li)
    vector=pd.DataFrame(vector)
    vector=np.transpose(vector)
    vektorNeg = vektorNeg.append(vector, ignore_index=True)
    #ilgili protein için k-mer dokümanı listeye eklemir
    kmer_listNeg.append(prot_k_mer)

vektorNeg=pd.DataFrame(vektorNeg)
writer = pd.ExcelWriter(r'C:\Users\zeyne\Desktop\Negatif.xlsx')
# write dataframe to excel
vektorNeg.to_excel(writer)
writer.save()
    
    
#EĞİTİLEN DOC2VEC MODELİ İLE TEST PROTEİNLERİN VEKTÖRLERE DÖNÜŞTÜRÜLMESİ
testVeri = pd.read_csv(r'C:\Users\zeyne\Desktop\Covid-PPI Proje\Veriler\HomoSapProteins.csv')
testVeri=testVeri['Entry']#sekans sütununun başlığı dosyada 'Entry'
print (df['Entry'])

k=3#k aminoasitlik birimler birer kelime olacak
sekansSayi=len(testVeri)
kmer_listtestVeri = []#k-merler ile ifade edilmiş pozitif protein dokümanlarını tutacak bu liste
vektortestVeri=[]
vektortestVeri=pd.DataFrame(vektortestVeri)
for i in range(sekansSayi):
    protein=testVeri[i]
    start=0 # start ve end değişkenleri substringin başlangıç ve bitiş indisleri. 
    end=start+k
    prot_k_mer="" #sıradaki proteinin k-merlerine bölünmüş halini bu değişken tutar. 
    k_merSayi=len(protein)-k+1
    j=0
    while (j<k_merSayi):
        k_mer=protein[start:end:1]
        prot_k_mer=prot_k_mer+" "+k_mer
        start=start+1
        end=start+k
        j=j+1
    li = list(prot_k_mer.split(" "))
    vector = modelDM3mer.infer_vector(li)
    vector=pd.DataFrame(vector)
    vector=np.transpose(vector)
    vektortestVeri = vektortestVeri.append(vector, ignore_index=True)
    #ilgili protein için k-mer dokümanı listeye eklemir
    kmer_listtestVeri.append(prot_k_mer)
    
vektortestVeri=pd.DataFrame(vektortestVeri)
vektortestVeri1=pd.DataFrame(vektortestVeri.iloc[0:100000,:])
vektortestVeri2=pd.DataFrame(vektortestVeri.iloc[100001:204185,:])

writer = pd.ExcelWriter(r'C:\Users\zeyne\Desktop\negData3Doc2Vec3Mer2.xlsx')
#write dataframe to excel
vektortestVeri2.to_excel(writer)
writer.save()

#NegTest=vektorNeg.iloc[NegRowNumTest,:] 
#NegTrain=vektorNeg.iloc[NegRowNumTrain,:]

now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)



#kod ile benzer doküman bulma
#modelBow3mer.docvecs.most_similar(positive=[modelBow3mer.infer_vector(li)],topn=5)

from sklearn.metrics.pairwise import cosine_similarity
#Gerçekte pozitif olan test kümesinin pozitif eğitim kümesindeki örneklere cos similarity ölçülür.
#negatif olmayan örneklerin ortalaması alınır
cos1 = cosine_similarity(PozTrain, PozTest)
ort1=[]
for i in range (len(PozTest)):
    toplam=0
    sayi=0
    for j in range(len(PozTrain)):
        if (cos1[j,i]>0):
            toplam=toplam+cos1[j,i]
            sayi=sayi+1
    ort1.append(toplam/sayi)
    
#Gerçekte pozitif olan test kümesinin negatşf eğitim kümesindeki örneklere cos similarity ölçülür.
#negatif olmayan örneklerin ortalaması alınır
cos2 = cosine_similarity(NegTrain, PozTest)
ort2=[]
for i in range (len(PozTest)):
    toplam=0
    sayi=0
    for j in range(len(NegTrain)):
        if (cos2[j,i]>0):
            toplam=toplam+cos2[j,i]
            sayi=sayi+1
    ort2.append(toplam/sayi)

#ortalama vektörleri posizitf test kümesindeki herbir örnek için sırasıyla pozitif ve negatif
#eğitim kümesindeki örneklere cos benzerliklerinin ortalamasını içerir. i. test örneği için
#ort1 değeri ort2 değerinden büyükse i. test örneğinin pozitif örneklere benzerliği daha fazladır
#Etiketi 1 olarak atanır. Aksi halde de 0 olacaktır.
tahmin1=[]
for i in range(len(PozTest)):
    if (ort1[i]>ort2[i]):
        tahmin1.append(1)
    else:
        tahmin1.append(0)

#aynı işlem negatif test örneklerine uygulanır.
cos1 = cosine_similarity(PozTrain, NegTest)
ort1=[]
for i in range (len(NegTest)):
    toplam=0
    sayi=0
    for j in range(len(PozTrain)):
        if (cos1[j,i]>0):
            toplam=toplam+cos1[j,i]
            sayi=sayi+1
    ort1.append(toplam/sayi)
    

cos2 = cosine_similarity(NegTrain, NegTest)
ort2=[]
for i in range (len(NegTest)):
    toplam=0
    sayi=0
    for j in range(len(PozTrain)):
        if (cos2[j,i]>0):
            toplam=toplam+cos2[j,i]
            sayi=sayi+1
    ort2.append(toplam/sayi)

tahmin2=[]
for i in range(len(NegTest)):
    if (ort1[i]>ort2[i]):
        tahmin2.append(1)
    else:
        tahmin2.append(0)

#tahmin1 gerçek sınıfı pozitif olan test örneklerinin tahmin edilen sınıflarıdır. tahmin 2 ise 
#negatif test örneklerinindir. tahmin1'de 1 değeri, tahmin2'de 0 değeri ne kadar fazlaysa o kadar
#iyidir. Buna göre accuracy hesaplanır.
x=tahmin1.count(1)
y=tahmin2.count(0)
acc=(x+y)/(67*2)

#Eğitim ve kümesinin etiketleri oluşturuldu. Proteinlerin olduğu excele eğitim ve test kümesindeki
#örneklere göre birer sütun vektörü tanımlayıp ordan okuma yaptım
trainLabel= pd.read_excel(r'C:\Users\zeyne\Desktop\Kaggle Data\Proje2\Veriler\Proteinler.xlsx', sheet_name='Sayfa3', header=0)
trainLabel=trainLabel['trLabel']
testLabel= pd.read_excel(r'C:\Users\zeyne\Desktop\Kaggle Data\Proje2\Veriler\Proteinler.xlsx', sheet_name='Sayfa4', header=0)
testLabel=testLabel['testLabel']

#eğitim ve test kümeleri oluşturuldu
train = pd.concat([PozTrain, NegTrain], axis=0)
test = pd.concat([PozTest, NegTest], axis=0)
           



model=models.Sequential()
model.add(layers.LSTM(28,input_shape=(X.shape[0], X.shape[1]),return_sequences=True))
model.add(layers.Dropout(0.4))
model.add(layers.LSTM(14))
model.add(layers.Dropout(0.5))
model.add(layers.Dense(2,activation="sigmoid"))
model.compile(optimizer='Adam',loss='binary_crossentropy',metrics=['accuracy'])
model.summary()
model.fit(X, Y,epochs=20,batch_size=512)
model.evaluate(test, testLabel)
len(Y)
Y.shape


vertical_stack = pd.concat([vektorPoz, vektorNeg], axis=0)


writer = pd.ExcelWriter(r'C:\Users\zeyne\Desktop\Mer.xlsx')
# write dataframe to excel
vertical_stack.to_excel(writer)
writer.save()


