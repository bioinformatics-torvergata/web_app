#!/usr/bin/env python
# coding: utf-8

# #Analisi di differenziale di espressione di geni/mirna/proteine
# #Gene Expression Quantification HTSeq - FPKM (17.753 samples)
# 

# Funzione che:
# - dato un gene/miRNA/protein
# - dato il nome del tumore (es. BRCA)
# - data una feature
# 
# Ci restituisce due boxplot per confrontare i due gruppi di pazienti con quella determinata feature.
# Calcoli il p-value con il test di Wilcoxon.
# 

 #Analisi di differenziale di espressione di geni/mirna/proteine
# #Gene Expression Quantification HTSeq - FPKM (17.753 samples)
# 

# Funzione che:
# - dato un gene/miRNA/protein
# - dato il nome del tumore (es. BRCA)
# - data una feature
# 
# Ci restituisce due boxplot per confrontare i due gruppi di pazienti con quella determinata feature.
# Calcoli il p-value con il test di Wilcoxon.
# 

import pandas as pd
import os 
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ranksums
import sys
import subprocess

import plotly.express as px
from function_new import open_dataframe_gene,df_feature,crealista
os.chdir("/Users/chiaranotturnogranieri/Downloads/notturno") #aggiungere in un file di configurazione questo path


if __name__ == "__main__":
    #parametri da inserire:
    gene= sys.argv[1]
    tumor= sys.argv[2]
    feature= sys.argv[3]

    cartella=sys.argv[4]
    #cartella='results'

    #df di interesse:
    ogg=open_dataframe_gene(gene,tumor,feature,cartella)

    if ogg!=0:
        df=ogg[0]
        dffeat=df_feature(ogg[1], tumor, feature) #df features
        
        if feature== 'patient_status':
            listaf0= crealista(dffeat, df,feature)
            
            lista_exp=[]#lista valori espressione
            lista_feature=[]
            for ele in df[listaf0]:
                lista_exp.append(df.loc[gene,ele])
                if int(ele[-3:-1])<11:
                    lista_feature.append('Tumor')
                else:
                    lista_feature.append('Ctrl')

            #controllo di avere sia dati tumorali che di controllo
            t=0
            c=0
            for ele in lista_feature:
                if ele== 'Tumor':
                    t+=1
                else:
                    c+=1
            if t==0 or c==0:
                print('non ci sono abbastanza dati')               
            
            
            
            d=pd.DataFrame({'nome01':listaf0, gene: lista_exp, feature: lista_feature})
        
        
        else:
            listaf, listaf01= crealista(dffeat, df,feature)
            dffeat=dffeat[['bcr_patient_barcode',feature]]

            #creo un dframe finale
            dframe=pd.DataFrame({'nome01':listaf01, 'bcr_patient_barcode':listaf})
            d= dframe.merge(dffeat)

            lista_exp=[] #lista valori espressione
            for ele in d['nome01']:
                lista_exp.append(df.loc[gene,ele])

            d[gene]=lista_exp

        #plt
        fig = px.scatter(x=range(10), y=range(10))
        fig=px.box(d,y=gene,x=feature,color=feature,points = 'all')
        if ogg[1]== "miRNA":
                fig.update_layout(yaxis_type="log")
        #fig.show()
      
        fig.write_html(cartella+'/'+gene+'_'+feature+'.html')
        
       

        #ranksum test per p-value
        p=list(set(d[feature]))
        
        df1_mask=d[feature]== p[0]
        dp0=d[df1_mask]
        
        df1_mask=d[feature]== p[1]
        dp1=d[df1_mask]

        w, p = ranksums(list(dp0[gene]), list(dp1[gene]))

        print(p)

