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
from function_new_prova import dataframe_OStime,overall_survival_analysis_pathway,open_gsva_df




if __name__ == "__main__":
    tumor= sys.argv[1]
    pathway= sys.argv[2]
    cartella=sys.argv[3]
    column=sys.argv[4]


    #open df con PASs(pathway activity scores)
    df=open_gsva_df(tumor)
  

    ###df con dati OS time
    OS=dataframe_OStime(tumor,column)



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

    overall_survival_analysis_pathway(pathway,tumor,df1,OS1,cartella)


