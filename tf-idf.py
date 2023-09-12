# -*- coding: utf-8 -*-
"""
Created on Sat Apr  9 15:03:17 2022

@author: zeynep
"""
import pandas as pd
import numpy as np

#UPLOAD POSITIVE PROTEINS
dfPozitif = pd.read_excel(r'Proteinler.xlsx', sheet_name='Human Proteins')
PozProt=dfPozitif['Protein Sequence']

#UPLOAD NEGATIVE PROTEINS
dfNegatif = pd.read_excel(r'negatifData.xlsx', sheet_name='Sayfa1')
NegProt=dfNegatif['sekans']

from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

#SPLIT POSITIVE PROTEINS TO THEIR K-MERS
k=3#Should repeat for different k value
sekansSayi=len(PozProt)
kmer_list = []
for i in range(sekansSayi):
    protein=PozProt[i]
    start=0  
    end=start+k
    prot_k_mer="" 
    k_merSayi=len(protein)-k+1
    j=0
    while (j<k_merSayi):
        k_mer=protein[start:end:1]
        prot_k_mer=prot_k_mer+" "+k_mer
        #prot_k_mer=prot_k_mer+" "+"'"+k_mer+"'"
        start=start+1
        end=start+k
        j=j+1
    #vektorPoz[i,1] = prot_k_mer
    kmer_list.append(prot_k_mer)
len(prot_k_mer)

#SPLIT NEGATIVE PROTEINS TO THEIR K-MERS
k=3#Should repeat for different k value
sekansSayi=len(NegProt)
#kmer_listNeg = []
for i in range(sekansSayi):
    protein=NegProt[i]
    start=0  
    end=start+k
    prot_k_mer=""  
    k_merSayi=len(protein)-k+1
    j=0
    while (j<k_merSayi):
        k_mer=protein[start:end:1]
        prot_k_mer=prot_k_mer+" "+k_mer
        start=start+1
        end=start+k
        j=j+1
    kmer_list.append(prot_k_mer)
len(kmer_list)  
len(prot_k_mer)
    

#CREATE VOCABULARY
#EVERY UNIQUE K-MER IS AN ENTRY OF THE VOCABULARY
i=0
protein=kmer_list[i]
sozluk=protein.split(" ")
for i in range (1,len(kmer_list)):
    protein=kmer_list[i]
    protein=protein.lstrip()
    kelime=protein.split(" ")
    sozluk= set(sozluk).union(set(kelime))
len(sozluk)

#CREATE COUNT MATRIX. 
countMatris=pd.DataFrame()
#THE NUMBER OF PASSES FOR POSITIVE AND NEGATIVE PROTEINS  IN THE DICTIONARY IS CALCULATED. 
for j in range (len(kmer_list)):
    protein=kmer_list[j]
    protein=protein.lstrip()
    kelime=protein.split(" ")
    wordDict = dict.fromkeys(sozluk, 0) 
    for i in range(len(kelime)):
        word=kelime[i]
        wordDict[word]+=1
    countMatris=countMatris.append(wordDict, ignore_index=True) 



countMatris.head
countMatris=countMatris.iloc[:,1:len(sozluk)]

from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

#CREATING TF-IDF MATRIX
import math #for logarithm function
tfMatris=pd.DataFrame()
j=0
i=0
N=len(kmer_list)
for i in range(1,N):
    a=countMatris.iloc[i,]
    tfVektor=np.zeros(len(sozluk))
    for j in range(len(sozluk)-1):
        tf=a[j]/(sum(a))
        df=np.count_nonzero(countMatris.iloc[:,j])
        idf=math.log((N/df),10)
        tfVektor[j]=tf*idf
    tfVektor=pd.DataFrame(tfVektor)
    tfVektor=np.transpose(tfVektor)
    tfMatris=tfMatris.append(tfVektor)

from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)


Label= pd.read_excel(rProteinler.xlsx', sheet_name='Sayfa3', header=0)
Label=Label['trLabel']

data2=tfMatris.insert(0, "etiket", Label)
data2.head
data2=pd.DataFrame(tfMatris)
writer = pd.ExcelWriter(r'tfMatrisEgitim.xlsx')
tfMatris.to_excel(writer)
writer.save()




#same process repeated for test dataset





