#!/usr/bin/env python
# coding: utf-8

# Differential analisys all tumor boxplot

# Passare un gene/miRNA/proteina e una feature


import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import ranksums
import subprocess
import os
import sys
import configparser
from pathlib import Path
from function_new_prova import read_clinical_data,open_dataframe_gene_boxplot_all_tumor,box_plot_all_tumor,df_feature_age,crealista_patient_status,pulisci_df,p_value,gene_mirna_proteina



#################################    
#             MAIN              #
#################################   

if __name__ == "__main__":
    #parametri in input:
    gene= sys.argv[1]
    feature= sys.argv[2]
    cartella=sys.argv[3]
    control=sys.argv[4]
    x=read_clinical_data()
    
    gene_input=gene
    if feature=="age_at_initial_pathologic_diagnosis":
        x1=df_feature_age(x,feature)
    else:
        x1=x

    #aprire df
    ogg_analisi=gene_mirna_proteina(gene, cartella,control) #es. miRNA, miRNA_ID, 'pathdataframe.csv'
    
    if ogg_analisi!=0:
        gene=ogg_analisi[3]
        #lista samples df expression gene/miRNA/prot
        dd=pd.read_csv(ogg_analisi[2],nrows=1)
        listamiRNA=list(dd.columns[1:])


        #patient status
        if feature== 'patient_status':    
            listanomi01, listanomi= crealista_patient_status(x, listamiRNA,ogg_analisi)
            
            #aprire df espressione
            d=open_dataframe_gene_boxplot_all_tumor(ogg_analisi[0],listanomi01, ogg_analisi[2],ogg_analisi[1]) #piu index
            

            lista_exp=[]#lista valori espressione
            lista_feature=[]#lista valori feature
            for ele in d.columns:
                lista_exp.append(d.loc[gene,ele])
                if int(ele[-3:-1])<11:
                    lista_feature.append('Tumor')
                else:
                    lista_feature.append('Ctrl')
                    

        else:
            #list samples with clinical data
            listap=list(x1.index) 
            
            listanomi=[] #nomi con cui ricavare i dati clinici
            listanomi01=[ogg_analisi[1]] #nomi con cui ricavare dati di espressione
            for el in listamiRNA:
                if el[:-4] in listap and el[:-4] not in listanomi:
                    if int(el[-3:-1])<11:
                        listanomi.append(el[:-4])
                        listanomi01.append(el)


            #df espressione solo dei campioni di cui abbiamo i dati clinici
            d=open_dataframe_gene_boxplot_all_tumor(ogg_analisi[0],listanomi01, ogg_analisi[2],ogg_analisi[1])



        colonna1=list(d.loc[gene,])
        colonna2=list(x1.loc[listanomi, "acronym"])
        if feature== 'patient_status':
            colonna3=lista_feature
        else:
            colonna3=list(x1.loc[listanomi, feature])


        df=pd.DataFrame({gene:colonna1, 'tumor':colonna2, feature:colonna3})
        df=df.dropna()


        #rimuovere i tumori che non hanno il ctrl 
        if feature=='patient_status':
            df=pulisci_df(df,feature)  



        box_plot_all_tumor(df, cartella, gene, feature,ogg_analisi[0],gene_input)
        p_value(df, cartella,feature,gene)

    else:
        print(0)
