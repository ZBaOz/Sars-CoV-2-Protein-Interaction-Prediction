# -*- coding: utf-8 -*-
"""
Created on Sat Apr  9 15:03:17 2022

@author: zeyne
"""
import pandas as pd
import numpy as np

#POZİTİF PROTEİNLERİN YÜKLENMESİ
dfPozitif = pd.read_excel(r'C:\Users\zeyne\Desktop\Covid-PPI Proje\Veriler\Proteinler.xlsx', sheet_name='Human Proteins')
PozProt=dfPozitif['Protein Sequence']

#NEGATİF PROTEİNLERİN YÜKLENMESİ
dfNegatif = pd.read_excel(r'C:\Users\zeyne\Desktop\Covid-PPI Proje\Veriler\negatifData3.xlsx', sheet_name='Sayfa1')
NegProt=dfNegatif['sekans']

from datetime import datetime
#çalışma süresi hesaplamak için get currenttime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

#POZİTİF PROTEİNLER K-MERLERİNE AYRILIR
k=3#k aminoasitlik birimler birer kelime olacak
sekansSayi=len(PozProt)
kmer_list = []#k-merler ile ifade edilmiş pozitif protein dokümanlarını tutacak bu liste
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
        #prot_k_mer=prot_k_mer+" "+"'"+k_mer+"'"
        start=start+1
        end=start+k
        j=j+1
    #ilgili protein için k-mer dokümanı listeye eklemir
    #vektorPoz[i,1] = prot_k_mer
    kmer_list.append(prot_k_mer)
len(prot_k_mer)

#NEGATİF PROTEİNLER K-MERLERİNE AYRILIR.
k=3#k aminoasitlik birimler birer kelime olacak
sekansSayi=len(NegProt)
#kmer_listNeg = []#k-merler ile ifade edilmiş negatif protein dokümanlarını tutacak bu liste
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
    #ilgili protein için k-mer dokümanı listeye eklemir
    kmer_list.append(prot_k_mer)
len(kmer_list)  
len(prot_k_mer)
    
#pozitif ve negatif proteinler tek bir k_mer listesinde:kmer_list

#SÖZLÜK OLUŞTURMA
#HER K-MERİNE AYRILMIŞ PROTEİN, 'split' İLE BOŞLUK KARAKTERİNE GÖRE KELİMELERE AYRILIR. 'SÖZLÜK' DEĞİŞKENİNE EKLENİR.
#EKLEME SIRASINDA VAR OLAN BİR K-MERİ YENİDEN EKLEMESİN DİYE 'union' FONKSİYONU KULLANILDI
i=0
protein=kmer_list[i]
sozluk=protein.split(" ")
for i in range (1,len(kmer_list)):
    protein=kmer_list[i]
    protein=protein.lstrip()#stringin başındaki boşluğu siler
    kelime=protein.split(" ")
    sozluk= set(sozluk).union(set(kelime))
#sözlük boyutu
len(sozluk)

#COUNT MATRİS OLUŞTURMA
countMatris=pd.DataFrame()
#POZİTİF VE NEGATİF PROTEİNLER İÇİN SÖZLÜKTEKİ KELİMELERİN GEÇME SAYILARI HESAPLANIR. 'countMatris'e EKLENİR.
for j in range (len(kmer_list)):
    protein=kmer_list[j]
    protein=protein.lstrip()#stringin başındaki boşluğu siler
    kelime=protein.split(" ")
    #wordDict sıradaki protein için frekans matrisini tutacaktır.Sözlük boyutu kadar boyuttan oluşur ve tüm değerler 0
    #ile ilklendirilir.
    wordDict = dict.fromkeys(sozluk, 0) 
    for i in range(len(kelime)):#bir kelime bir proteinde geçtikçe ilgili sütunun değeri 1 artırılır
        word=kelime[i]
        wordDict[word]+=1
    countMatris=countMatris.append(wordDict, ignore_index=True) 



countMatris.head
#sözlük oluşturuken boşluğu da bir karakter olarak bir defa alıyor ve bu count matrisin ilk sütunu. tf idf oluştururken 
#0'a bölme hatası oluşturuyor bu yüzden sildim
countMatris=countMatris.iloc[:,1:len(sozluk)]

from datetime import datetime
#çalışma süresi hesaplamak için get currenttime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

#TF-IDF MATRİSİ OLUŞTURMA
import math #logaritma fonksiyonu için
tfMatris=pd.DataFrame()#tf-idf değerlerini tutacak olan matris
j=0
i=0
N=len(kmer_list)#toplam doküman/protein sayısı
for i in range(1,N):
    a=countMatris.iloc[i,]#count matristeki i. satır vektörü a değişkenine alınır. 
    tfVektor=np.zeros(len(sozluk)) #tfVektor i. protein için tf-idf degerlerini içerecek satır vektörü
    for j in range(len(sozluk)-1):#a vektörünün tüm elemanlarını dolaşacak. Tüm proteinler için boyut sözlük boyutu kadardır
        #tf=sıradaki kelimenin i. proteinde geçme sayısı/proteindeki toplam kelime sayısı. a vektöründe her bir kelimenin
        #o proteinde geçme sayısı olduğundan bu sayıların toplamı proteindeki kelime sayısını verir.
        tf=a[j]/(sum(a))
        #df=j. kelimenin tüm dokümanlarda geçme sayısıdır.Bir dokümanda 1 defa geçmesi ile 5 defa geçmesi arasında bir
        #fark yok. dokümanda geçip geçmediği önemli. bu nedenle 0 olmayan hücre sayısı hesaplanır. Bu kaç proteinde geçtiğini
        #verecektir.
        df=np.count_nonzero(countMatris.iloc[:,j])
        #idf=toplam doküman sayısı/ilgili kelimenin geçtiği doküman sayısının 10 luk tabanda logaritmasıdır
        idf=math.log((N/df),10)
        tfVektor[j]=tf*idf
    #tf vektörü i. proteinin tf-idf değerlerini içerir. vektör data frame e çevrildiğinde sütun vektörü olduğundan
    #transpozu alınarak tfMatrisine eklenir.
    tfVektor=pd.DataFrame(tfVektor)
    tfVektor=np.transpose(tfVektor)
    tfMatris=tfMatris.append(tfVektor)

from datetime import datetime
#çalışma süresi hesaplamak için get currenttime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

#Eğitim ve kümesinin etiketleri oluşturuldu. Proteinlerin olduğu excele eğitim ve test kümesindeki
#örneklere göre birer sütun vektörü tanımlayıp ordan okuma yaptım
Label= pd.read_excel(r'C:\Users\zeyne\Desktop\Covid-PPI Proje\Veriler\Proteinler.xlsx', sheet_name='Sayfa3', header=0)
Label=Label['trLabel']

#tfmatrisine etiket bilgisi eklenir ve kaydedilir.
data2=tfMatris.insert(0, "etiket", Label)
data2.head
data2=pd.DataFrame(tfMatris)
#excel yazdırma
writer = pd.ExcelWriter(r'C:\Users\zeyne\Desktop\tfMatrisEgitim.xlsx')
# write dataframe to excel
tfMatris.to_excel(writer)
writer.save()


#4 mer için tf matrisini matlabde oluşturdum. Burada uzun sürüyor diye. Oyüzden sadece count matrisi kaydettim
countTest=pd.DataFrame(countMatrisTest)
countTest.to_csv(r'C:\Users\zeyne\Desktop\countMatrisTest.csv')

#kMerList yazdırılır
data3=pd.DataFrame(kmer_list)   
writer = pd.ExcelWriter(r'C:\Users\zeyne\Desktop\kMerListEgitim.xlsx')
# write dataframe to excel
data3.to_excel(writer)
writer.save() 


#test verisinin aynı formatta oluşturulması için aday proteinler yüklenir
df = pd.read_csv(r'C:\Users\zeyne\Desktop\Covid-PPI Proje\Veriler\HomoSapSekans.csv')
sekanslar=df['V1']#sekans sütununun başlığı dosyada 'Entry'
print (df['V1'])

#TEST VERİSİ K-MERLERİNE AYRILIR
k=3#k aminoasitlik birimler birer kelime olacak
sekansSayi=len(sekanslar)
kmer_listTest = []#k-merler ile ifade edilmiş pozitif protein dokümanlarını tutacak bu liste
for i in range(sekansSayi):
    protein=sekanslar[i]
    start=0 # start ve end değişkenleri substringin başlangıç ve bitiş indisleri. 
    end=start+k
    prot_k_mer="" #sıradaki proteinin k-merlerine bölünmüş halini bu değişken tutar. 
    k_merSayi=len(protein)-k+1
    j=0
    while (j<k_merSayi):
        k_mer=protein[start:end:1]
        prot_k_mer=prot_k_mer+" "+k_mer
        #prot_k_mer=prot_k_mer+" "+"'"+k_mer+"'"
        start=start+1
        end=start+k
        j=j+1
    #ilgili protein için k-mer dokümanı listeye eklemir
    #vektorPoz[i,1] = prot_k_mer
    kmer_listTest.append(prot_k_mer)
len(kmer_listTest)




#TEST VERİSİNİN COUNT MATRİSİ OLUŞTURULUR
countMatrisTest=pd.DataFrame()
#TEST KÜMESİNDEKİ PROTEİNLER İÇİN SÖZLÜKTEKİ KELİMELERİN GEÇME SAYILARI HESAPLANIR. 'countMatris'e EKLENİR.
#for j in range (len(kmer_listTest)):
for j in range (132507,len(kmer_listTest)):
    protein=kmer_listTest[j]
    protein=protein.lstrip()#stringin başındaki boşluğu siler
    kelime=protein.split(" ")
    #wordDict sıradaki protein için frekans matrisini tutacaktır.Sözlük boyutu kadar boyuttan oluşur ve tüm değerler 0
    #ile ilklendirilir.
    wordDict = dict.fromkeys(sozluk, 0) 
    for i in range(len(kelime)):#bir kelime bir proteinde geçtikçe ilgili sütunun değeri 1 artırılır
        word=kelime[i]
        if word in wordDict:
            wordDict[word]+=1
    countMatrisTest=countMatrisTest.append(wordDict, ignore_index=True) 



wordDict[word].

count1=pd.DataFrame(countMatris.iloc[:,1:70000])

count1=np.transpose(countMatris)

count2=count1.iloc[20000:40000,:]
count2.head

sozluk2=pd.DataFrame(sozluk)
writer = pd.ExcelWriter(r'C:\Users\zeyne\Desktop\sozluk.xlsx')
# write dataframe to excel
sozluk2.to_excel(writer)
writer.save()

import csv
with open(r'C:\Users\zeyne\Desktop\countMatris.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerows(countMatris)

countMatris.to_csv(r'C:\Users\zeyne\Desktop\countMatris.csv')

df2 = pd.DataFrame(X.toarray(), columns=tf.get_feature_names())
print(df2)
df2.to_csv('C:/Users/zeyne/Desktop/data.csv', index=False)
df2.head

with open ("tf_idf.csv",'a', newline='') as file:
    writer = csv.writer(file)
    writer.writerow([df2])
    file.close()






