#!/usr/bin/env python
# coding: utf-8

# ANALISI OVERALL SURVIVAL
# 


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import sys
import os
from function_new import open_dataframe_gene_overall,dataframe_OStime,overall_survival_analysis



if __name__ == "__main__":
    #dfclinic=read_clinical_data()
    #dfclinic=pd.read_csv('/mnt/data/notturno/TCGA-CDR-SupplementalTableS1.csv')
    gene= sys.argv[1] #gene/mirna/protein
    tumor=sys.argv[2]
    cartella=sys.argv[3]
    feature='median'



    df=open_dataframe_gene_overall(gene,tumor)

    OS=dataframe_OStime(tumor)

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

    overall_survival_analysis(gene,tumor,feature, cartella,df1,OS1,gene)



