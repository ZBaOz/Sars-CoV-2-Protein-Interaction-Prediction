# -*- coding: utf-8 -*-
"""
Created on Fri Apr 15 20:25:47 2022

@author: zeynep
"""
import pandas as pd
import numpy as np

pip install xlrd
import xlrd
pip install openpyxl

#for training word2vec model, 204,185 human proteıns were downloaded from unıprot. 
df = pd.read_csv(r'HomoSapSekans.csv')
sekanslar=df['V1']
#print (df['Entry'])
sekanslar.head
%reset
#READ POSITIVE AND NEGATIVE HUMAN PROTEINS
tumData = pd.read_excel(r'VeriSeti3.xlsx', sheet_name='Sayfa1', engine='openpyxl')
PozProt=tumData['human sekans']
NegProt=tumData['negatif sekans']

#POSITIVE AND NEGATIVE HUMAN PROTEINS WERE ALSO ADDED TO TRAINING DATASET. THERE ARE 204,849 UNIQOE PROTEINS IN THE TRAINING DATASET.
tumData = pd.concat([sekanslar, PozProt], axis=0)
tumData = pd.concat([tumData, NegProt], axis=0)

tumData.head

#LIBRARY FOR TOKENIZATION
import nltk 
WPT = nltk.WordPunctTokenizer()

from datetime import datetime
now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

#PROTEINS IN TRAINING DATASET SPLITTED TO THEIRS K-MERS.
k=3#shoul repeat for different k values
sekansSayi=len(tumData)
kmer_list = []
for i in range(sekansSayi):
    protein=tumData.iloc[i]
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
    prot_k_mer=prot_k_mer.lstrip()
    tokens=WPT.tokenize(prot_k_mer)
    kmer_list.append(tokens)
len(kmer_list)    

from gensim.models import Word2Vec

#size=the size of the output vector
#workers=this parameter is for parallelization.
#window=Maximum distance between the current and predicted word within a sentence.
#min_count= Ignores all words with total frequency lower than this.
#sg ({0, 1}– Training algorithm: 1 for skip-gram; otherwise CBOW. 
#hs ({0, 1}, optional) – If 1, hierarchical softmax will be used for model training. If 0, and negative is non-zero, negative sampling will be used.
#cbow_mean ({0, 1}, optional) – If 0, use the sum of the context word vectors. If 1, use the mean, only applies when cbow is used.
model3Mer = Word2Vec(kmer_list, min_count=1, vector_size=300, window=10, sg=0, negative=5, epochs=20,hs=0)
model3Mer.save(r"C:\Users\zeyne\Desktop\model3Mer.model")

word_vectors = model3Mer.wv
word_vectors.save(r"word2vec.wordvectors")
#print(model3Mer)
  
#POSITIVE PROTEINS ALSO SPLITTED TO THEIRS K-MERS.
k=3
sekansSayi=len(PozProt)
Pozkmer_list = []
for i in range(sekansSayi):
    protein=PozProt.iloc[i]
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
    prot_k_mer=prot_k_mer.lstrip()
    tokens=WPT.tokenize(prot_k_mer)
    Pozkmer_list.append(tokens)
    
#len(Pozkmer_list)   

data=pd.DataFrame()#poz ve neg veriyi tutacak olan değişken

#Feature vectors for positive protins obatined using word2vec model. 
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

#NEGATIVE PROTEINS ALSO SPLITTED TO THEIRS K-MERS.
k=3
sekansSayi=len(NegProt)
Negkmer_list = []
for i in range(sekansSayi):
    protein=NegProt.iloc[i]
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
    prot_k_mer=prot_k_mer.lstrip()
    tokens=WPT.tokenize(prot_k_mer)
    Negkmer_list.append(tokens)


for i in range (len(Negkmer_list)):
    tokens=Negkmer_list[i]
    sent_list=[]
    for word in tokens:
        wv2=model3Mer.wv[word]
        sent_list.append(wv2)
    vektor=np.mean(sent_list,axis=0)
    vektor=pd.DataFrame(vektor)
    vektor=np.transpose(vektor)
    data = data.append(vektor, ignore_index=True)


#TEST PROTEINS ALSO SPLITTED TO THEIRS K-MERS.
test = pd.read_csv(r'HomoSapSekans.csv')
test=test['V1']'
print (test)
test.head
#test protein sekansları k-merlerine ayır
k=3
sekansSayi=len(test)
testkmer_list = []
for i in range(sekansSayi):
    protein=test.iloc[i]
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
    prot_k_mer=prot_k_mer.lstrip()
    tokens=WPT.tokenize(prot_k_mer)
    testkmer_list.append(tokens)

data=pd.DataFrame()


#for i in range (len(testkmer_list)):
    tokens=testkmer_list[i]
    sent_list=[]
    for word in tokens:
        wv2=model5Mer.wv[word]
        sent_list.append(wv2)
    vektor=np.mean(sent_list,axis=0)
    vektor=pd.DataFrame(vektor)
    vektor=np.transpose(vektor)
    data = data.append(vektor, ignore_index=True)





