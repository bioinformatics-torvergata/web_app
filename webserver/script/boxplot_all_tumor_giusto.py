#!/usr/bin/env python
# coding: utf-8

# Differential analisys all tumor boxplot

# Passare un gene/miRNA/proteina e una feature

# In[1]:


import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import subprocess
import os
import sys


# In[9]:



def gene_mirna_proteina(gene):
    if gene in open("/mnt/data/notturno/miRNA/namemirna.txt").read().split("\n"):
        return(['miRNA',"miRNA_ID",'/mnt/data/notturno/miRNA/DataFrameTCGA_miRNA.csv'])
    
    if gene in open("/mnt/data/notturno/gene_expression/ENSG.txt").read().split("\n"):
        listageni=open("/mnt/data/notturno/gene_expression/ENSG.txt").read().strip().split("\n")
        posizione="-e "+str(listageni.index(gene)+1)+"p"
        
        #creiamo un dataframe piu piccolo dove c'è solo la riga del gene che è stato selezionato
        path=cartella+'/'+gene+"_df.txt"
        out_file=open(path,"w")
        subprocess.call(["sed","-n", "-e 1p", posizione, "/mnt/data/notturno/gene_expression/DataFrameTCGA_gene.csv"],stdout=out_file)
        
        return(['gene','gene_id',path])
    
    
    if gene in open("/mnt/data/notturno/protein/namepeptide.csv").read().split("\n"):
        return(['protein', "peptide_target",'/mnt/data/notturno/protein/DataFrameTCGA_protein.csv'])
    
    else:
        return('per il nome inserito non è disponibile la ricerca')


# In[3]:


def open_dataframe_gene(gene,listanomi01, path_dataframe):
    if gene == 'miRNA':
        df=pd.read_csv(path_dataframe,usecols=listanomi01)
        df=df.set_index('miRNA_ID')
        return(df)
    if gene == 'gene':
        df=pd.read_csv(path_dataframe,usecols=listanomi01)
        df=df.set_index("gene_id")
        return (df)
    if gene =='protein':
        df=pd.read_csv(path_dataframe,usecols=listanomi01)
        df=df.set_index('peptide_target')
        return (df)
    else: 
        print("per il nome inserito non è disponibile la ricerca")
        return 0


# In[13]:


def box_plot(df1, cartella):
    os.mkdir(cartella)
    sns.set(rc={'figure.figsize':(25.7,8.27)})
    sns.set_style("white")

    my_order = df1.groupby(by=["tumor"])[gene].median().iloc[::-1].sort_values().index

    ax=sns.boxplot(x="tumor", y=gene, hue=feature, data=df1, palette="Set2", width=0.7, order=my_order)
    ax.set_yscale("log")

    plt.savefig(cartella+'/'+gene+'_'+feature+'.png')


# In[5]:


def pulizia_df(df):
    return(df.dropna())


# In[19]:



path="/mnt/data/notturno/"

gene=sys.argv[1]
feature=sys.argv[2]

cartella=path+'web_app/webserver/rolls/static/media/saveanalisi/boxplot_all_tumor/'+ sys.argv[3]

x=pd.read_csv("/mnt/data/notturno/clinical_PANCAN_patient_with_followup_modificato.csv")
x1=x.set_index("bcr_patient_barcode")


#aprire df
ogg_analisi=gene_mirna_proteina(gene)
print(ogg_analisi)

d=pd.read_csv(ogg_analisi[2],nrows=1)
listamiRNA=list(d.columns[1:])
#print('lista pazienti miRNA',len(listamiRNA))

listap=list(x['bcr_patient_barcode']) #lista pazienti con dati clinici


listanomi=[] #nomi con cui ricavare i dati clinici
listanomi01=[ogg_analisi[1]] #nomi con cui ricavare dati di espressione
for el in listamiRNA:
    if el[:-4] in listap:
        listanomi.append(el[:-4])
        listanomi01.append(el)
print(len(listanomi01))
        
#aprire il dataframe di interesse
d=open_dataframe_gene(ogg_analisi[0],listanomi01, ogg_analisi[2])


colonna1=list(d.loc[gene,])
colonna2=list(x1.loc[listanomi, "acronym"])
colonna3=list(x1.loc[listanomi, feature])


df=pd.DataFrame({gene:colonna1, 'tumor':colonna2, feature:colonna3})
df=pulizia_df(df)

box_plot(df, cartella, gene, feature)


# In[ ]:




