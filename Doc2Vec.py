# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 19:39:17 2022

@author: zeyne
"""

import pandas as pd
import numpy as np
#204,185 human proteins were downloaded from uniprot to create the training dataset
df = pd.read_csv(r'\HomoSapProteins.csv')
sekanslar=df['Entry']#the caption of sequence column is 'entry'
print (df['Entry'])
df.head
%reset
#UPLOAD POSITIVE PROTEINS
dfPozitif = pd.read_excel(r'\Proteinler.xlsx', sheet_name='Human Proteins')
PozProt=dfPozitif['Protein Sequence']

#UPLOAD NEGATTIVE PROTEINS
dfNegatif = pd.read_excel(r'\negatifData3.xlsx', sheet_name='Sayfa1')
NegProt=dfNegatif['sekans']

#POSITIVE AND NEGATIVE PROTEINS ARE ALSO ADDED TO DOC2VEC TRAINSET. AFTER ADDITION THERE ARE 204.849 UNIQUE HUMAN PROTEINS
sekanslar = pd.concat([sekanslar, PozProt], axis=0)
sekanslar = pd.concat([sekanslar, NegProt], axis=0)

#SPLIT DATASET AS TRAIN AND TEST RANDOMLY
import random
indis=np.arange(start=0, stop=332, step=1)
PozRowNumTest=random.sample(range(0, 332), 67)
PozRowNumTrain=np.delete(indis, PozRowNumTest)
NegRowNumTest=random.sample(range(0, 332), 67)
NegRowNumTrain=np.delete(indis, NegRowNumTest)


from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

#SPLIT PROTEIN SEQUENCES TO THEIR K-MERS. AFTER THAT, ALL PROTEINS ARE CONSISTED OF WORDS LENGTH OF K AMINO ACIDS.
k=3 #THIS PROCESS WAS REPEATED FOR DIFFERENT K-VALUES
sekansSayi=len(sekanslar)
kmer_list = []
for i in range(sekansSayi):
    protein=sekanslar.iloc[i]
    start=0  
    end=start+k
    prot_k_mer="["  
    k_merSayi=len(protein)-k+1
    j=0
    while (j<k_merSayi-1):
        k_mer=protein[start:end:1]
        prot_k_mer=prot_k_mer+"'"+k_mer+"',"+" "
        start=start+1
        end=start+k
        j=j+1
    k_mer=protein[start:end:1]
    prot_k_mer=prot_k_mer+"'"+k_mer+"'"
    prot_k_mer=prot_k_mer+"]"
    #ADD RELATED K-MER PROTEIN TO THE LIST
    kmer_list.append(prot_k_mer)

len(kmer_list)
sekanslar.iloc[204183]
 
#TRAINING DOC2VEC MODEL
from gensim.test.utils import common_texts
from gensim.models.doc2vec import Doc2Vec, TaggedDocument

del documents

documents = [TaggedDocument(doc, [i]) for i, doc in enumerate(kmer_list)]
len(documents)
#doc2Vec function parameters: 
#dm:training algorithm. If the value is 1, training algorithm is distributed memory. 
#vector_size= The size of output vector
#window=The maximum distance between the current and predicted word within a sentence.
#Ignores all words with total frequency lower than this.
modelDM3mer = Doc2Vec(documents, vector_size=300, window=10, min_count=1, workers=4, dm=1, epochs=20, negative=5)


now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

#SAVE THE MODEL
modelDM3mer.save("modelDM3mer.mod")
vector = model4mer.infer_vector(["ILR", "IGN"])
model = Doc2Vec.load("./doc2vecmodel.mod")


#OBTAINING POSITIVE PROTEIN VECTORS BY DOC2VEC MODEL. 
k=3#THIS PROCESS WAS REPEATED FOR DIFFERENT K-VALUES
sekansSayi=len(PozProt)
kmer_listPoz = []
vektorPoz=[]
vektorPoz=pd.DataFrame(vektorPoz)
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
        start=start+1
        end=start+k
        j=j+1
    li = list(prot_k_mer.split(" "))
    vector = modelDM3mer.infer_vector(li)
    vector=pd.DataFrame(vector)
    vector=np.transpose(vector)
    vektorPoz = vektorPoz.append(vector, ignore_index=True)
    kmer_listPoz.append(prot_k_mer)


#PozTest=vektorPoz.iloc[PozRowNumTest,:] 
#PozTrain=vektorPoz.iloc[PozRowNumTrain,:] 
#PozTest.shape

vektorPoz=pd.DataFrame(vektorPoz)
writer = pd.ExcelWriter(r'Pozitif.xlsx')
#write dataframe to excel
vektorPoz.to_excel(writer)
writer.save()
    

#OBTAINING POSITIVE PROTEIN VECTORS BY DOC2VEC MODEL.
k=3#THIS PROCESS WAS REPEATED FOR DIFFERENT K-VALUES
sekansSayi=len(NegProt)
kmer_listNeg = []
vektorNeg=[]
vektorNeg=pd.DataFrame(vektorNeg)
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
    li = list(prot_k_mer.split(" "))
    vector = modelDM3mer.infer_vector(li)
    vector=pd.DataFrame(vector)
    vector=np.transpose(vector)
    vektorNeg = vektorNeg.append(vector, ignore_index=True)
    kmer_listNeg.append(prot_k_mer)

vektorNeg=pd.DataFrame(vektorNeg)
writer = pd.ExcelWriter(r'Negatif.xlsx')
# write dataframe to excel
vektorNeg.to_excel(writer)
writer.save()
    
    
#OBTAINING POSITIVE PROTEIN VECTORS BY DOC2VEC MODEL.
testVeri = pd.read_csv(r'\HomoSapProteins.csv')
testVeri=testVeri['Entry']
print (df['Entry'])

k=3
sekansSayi=len(testVeri)
kmer_listtestVeri = []
vektortestVeri=[]
vektortestVeri=pd.DataFrame(vektortestVeri)
for i in range(sekansSayi):
    protein=testVeri[i]
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
    li = list(prot_k_mer.split(" "))
    vector = modelDM3mer.infer_vector(li)
    vector=pd.DataFrame(vector)
    vector=np.transpose(vector)
    vektortestVeri = vektortestVeri.append(vector, ignore_index=True)
    kmer_listtestVeri.append(prot_k_mer)
    
vektortestVeri=pd.DataFrame(vektortestVeri)

writer = pd.ExcelWriter(r'negData3Doc2Vec3Mer2.xlsx')
#write dataframe to excel
vektortestVeri.to_excel(writer)
writer.save()


now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)




from sklearn.metrics.pairwise import cosine_similarity
#The cos similarity of the actually positive test set to the samples in the positive training set is measured.
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
    
#The cos similarity of the actually positive test set to the samples in the negative training set is measured.
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

#The mean vectors contain the average of the cos similarities for each sample in the positive test set to the samples in the positive and negative training set, respectively. 
#If the mean1 value for test sample i is greater than the mean2 value, test sample i is more similar to positive samples. Its label is assigned as 1. Otherwise it will be 0.
tahmin1=[]
for i in range(len(PozTest)):
    if (ort1[i]>ort2[i]):
        tahmin1.append(1)
    else:
        tahmin1.append(0)

#Same process is applied to negative proteins
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

#tahmin1 are the predicted classes of test samples whose true class is positive. tahmin 2 is for negative test samples. 
#The more values of 1 in tahmin1 and 0 in tahmin2, the better. Accuracy is calculated accordingly.
x=tahmin1.count(1)
y=tahmin2.count(0)
acc=(x+y)/(67*2)



