#!/usr/bin/env python
# coding: utf-8

# OS PATHWAY ACTIVITY SCORE
# 
# bisogna dare in input: 
# - tumore
# - nome del pathway
# 


import pandas as pd
import sys
import os
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr as prs
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines.plotting import add_at_risk_counts
from statsmodels.stats import multitest as multi
import numpy as np
from function_new_prova import dataframe_OStime


# def dataframe_OStime(dfclinic):
#     OS=(dfclinic[['bcr_patient_barcode','OS.time']]) 
#     df1_mask=dfclinic['type']==tumor
#     OS=dfclinic[df1_mask]
#     OS=(OS[['bcr_patient_barcode','OS.time']]) 
#     OS=OS.set_index('bcr_patient_barcode')
#     #display(OS)
#     OS=OS.dropna()
#     return(OS)



def overall_survival_analysis(m,tumor,df1,OS1,cartella):
      
    i1=df1.loc[m,:] > df1.loc[m,:].median()
    i2 = df1.loc[m,:] < df1.loc[m,:].median() 
    
    kmf = KaplanMeierFitter()
    

    
    results = logrank_test((OS1[i1]), (OS1[i2]),list(df1.loc[m,i1]),list(df1.loc[m,i2]), alpha=.95)
    
    if results.p_value < 1:
        os.mkdir(cartella)
        print("p-value:",results.p_value)
        
        kmf.fit((OS1[i1]), label="Higher expression")
        a1 = kmf.plot()

        kmf.fit((OS1[i2]), label="Lower expression")
        kmf.plot(ax=a1)
        
        plt.savefig(cartella+"/"+m+"_"+tumor+".png")

    else:
        print("pvalue>1")
   


tumor= sys.argv[1]
pathway= sys.argv[2]
cartella=sys.argv[3]
#cartella="/mnt/data/notturno/web_app/webserver/rolls/static/media/saveanalisi/"+ sys.argv[3]

#open df con PASs(pathway activity scores)
df=pd.read_csv('/mnt/data/notturno/gsva/gsva_'+tumor+'.csv')
df=df.set_index('Unnamed: 0')
df.columns= [x.replace(".","-") for x in df.columns]



###df con dati OS time
#dfclinic=pd.read_csv('/mnt/data/notturno/TCGA-CDR-SupplementalTableS1.csv')
OS=dataframe_OStime(tumor)



lista=list(OS.index)
oslist=[]
dflist=[]
for name in df.columns:
    if name[:-4] in lista:
        if name[-1]!="x" and name[-1]!="y":
            if int(name[-3:-1])<10:
                oslist.append(name[:-4]) #lista dei campioni del tumore di cui abbiamo OS.time
                dflist.append(name) #lista campioni 01

df1=df[dflist]


df1.columns= [(x[:-4]) for x in df1.columns]


OS1=OS.loc[oslist,:]

df1=df1*10

overall_survival_analysis(pathway,tumor,df1,OS1,cartella)

