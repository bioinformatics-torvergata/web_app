#!/usr/bin/env python
# coding: utf-8

# ANALISI OVERALL SURVIVAL
# 

# In[7]:


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr as prs
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines.plotting import add_at_risk_counts
from statsmodels.stats import multitest as multi
import numpy as np
import sys
import os
#from function_new import open_dataframe_gene_boxplot_all_tumor,read_clinical_data


def open_dataframe_gene(gene,tumor):
    if gene in open("/mnt/data/notturno/miRNA/namemirna.txt").read().split("\n"):
        df=pd.read_csv('/mnt/data/notturno/miRNA/DataFrameTCGA_miRNA.csv')
        df=df.set_index('miRNA_ID')
        return(df)
    if gene in open("/mnt/data/notturno/gene_expression/ENSG.txt").read().split("\n"):
        df=pd.read_csv("/mnt/data/notturno/Dataframe_tumorgene/Dataframe_FPKM/Dataframe_FPKM_"+tumor+".csv")
        df=df.set_index("gene_id")
        return (df)
    if gene in open("/mnt/data/notturno/protein/namepeptide.csv").read().split("\n"):
        df=pd.read_csv("/mnt/data/notturno/protein/DataFrameTCGA_protein.csv")
        df=df.set_index('peptide_target')
        return (df)
    else: 
        print("per il nome inserito non è disponibile la ricerca")
        return 0


# In[9]:


def dataframe_OStime(dfclinic):
    OS=(dfclinic[['bcr_patient_barcode','OS.time']]) 
    df1_mask=dfclinic['type']==tumor
    OS=dfclinic[df1_mask]
    OS=(OS[['bcr_patient_barcode','OS.time']]) 
    OS=OS.set_index('bcr_patient_barcode')
    #display(OS)
    OS=OS.dropna()
    return(OS)


# In[12]:


def overall_survival_analysis(m,tumor,feature,cartella,df1,OS1):
    
    
    i1=df1.loc[m,:] > df1.loc[m,:].median()
    i2 = df1.loc[m,:] < df1.loc[m,:].median() 
    
    kmf = KaplanMeierFitter()
    

    #if np.mean(list(df1.loc[m,i2]))>0:
    results = logrank_test((OS1[i1]), (OS1[i2]),list(df1.loc[m,i1]),list(df1.loc[m,i2]), alpha=.95)
    
    if results.p_value < 1:
        os.mkdir("../rolls/static/media/saveanalisi/"+cartella)
        print("p-value:",results.p_value)
        
        kmf.fit((OS1[i1]), list(df1.loc[m,i1]), label="Higher expression")
        a1 = kmf.plot()

        kmf.fit((OS1[i2]),list(df1.loc[m,i2]) , label="Lower expression")
        kmf.plot(ax=a1)
        plt.savefig("../rolls/static/media/saveanalisi/"+cartella+"/overallsurvival_"+gene+"_"+tumor+".png")
        
    else:
        print("ANALISI NON VALIDA pvalue>1")
    #else: 
     #   print("ANALISI NON VALIDA MEDIA <0 ??")


# In[ ]:





# In[13]:


dfclinic=pd.read_csv('/mnt/data/notturno/TCGA-CDR-SupplementalTableS1.csv')

#dfclinic=read_clinical_data()

gene= sys.argv[1] #gene/mirna/protein
tumor=sys.argv[2]
cartella=sys.argv[3]
feature='median'


df=open_dataframe_gene(gene,tumor)
OS=dataframe_OStime(dfclinic)

lista=list(OS.index)
oslist=[]
dflist=[]
for name in df:
     if name[:-4] in lista:
        if name[-1]!="x" and name[-1]!="y":
            if int(name[-3:-1])<10:
                oslist.append(name[:-4]) #lista dei campioni del tumore di cui abbiamo OS.time
                dflist.append(name)

df1=df[dflist]
df1.columns= [(x[:-4]) for x in df1.columns]


OS1=OS.loc[oslist,:]

overall_survival_analysis(gene,tumor,feature, cartella,df1,OS1)



