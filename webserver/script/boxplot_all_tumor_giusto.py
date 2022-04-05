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


get_ipython().run_cell_magic('time', '', 'def gene_mirna_proteina(gene):\n    if gene in open("/mnt/data/notturno/miRNA/namemirna.txt").read().split("\\n"):\n        return([\'miRNA\',"miRNA_ID",\'/mnt/data/notturno/miRNA/DataFrameTCGA_miRNA.csv\'])\n    \n    if gene in open("/mnt/data/notturno/gene_expression/ENSG.txt").read().split("\\n"):\n        listageni=open("/mnt/data/notturno/gene_expression/ENSG.txt").read().strip().split("\\n")\n        posizione="-e "+str(listageni.index(gene)+1)+"p"\n        \n        #creiamo un dataframe piu piccolo dove c\'è solo la riga del gene che è stato selezionato\n        path=cartella+gene+"_df.txt"\n        out_file=open(path,"w")\n        subprocess.call(["sed","-n", "-e 1p", posizione, "/mnt/data/notturno/gene_expression/DataFrameTCGA_gene.csv"],stdout=out_file)\n        \n        return([\'gene\',\'gene_id\',path])\n    \n    \n    if gene in open("/mnt/data/notturno/protein/namepeptide.csv").read().split("\\n"):\n        return([\'protein\', "peptide_target",\'/mnt/data/notturno/protein/DataFrameTCGA_protein.csv\'])\n    \n    else:\n        return(\'per il nome inserito non è disponibile la ricerca\')')


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


def box_plot(df1):
   
    sns.set(rc={'figure.figsize':(25.7,8.27)})
    sns.set_style("white")

    my_order = df1.groupby(by=["tumor"])[gene].median().iloc[::-1].sort_values().index

    ax=sns.boxplot(x="tumor", y=gene, hue=feature, data=df1, palette="Set2", width=0.7, order=my_order)
    ax.set_yscale("log")

    plt.show()
    #plt.savefig('save_as_a_png1log.png')


# In[5]:


def pulizia_df(df):
    return(df.dropna())


# In[19]:


get_ipython().run_cell_magic('time', '', 'path="/mnt/data/notturno/"\n#gene="ACVRL1"\n\n#gene="hsa-let-7a-1"\nfeature=\'menopause_status\'\n\ngene=\'ENSG00000167700.7\'\ncartella=\'provabox/\'            #sys.argv[3]\n\nx=pd.read_csv("/mnt/data/notturno/clinical_PANCAN_patient_with_followup_modificato.csv")\nx1=x.set_index("bcr_patient_barcode")\n\n\n#aprire df\nogg_analisi=gene_mirna_proteina(gene)\nprint(ogg_analisi)\n\nd=pd.read_csv(ogg_analisi[2],nrows=1)\nlistamiRNA=list(d.columns[1:])\n#print(\'lista pazienti miRNA\',len(listamiRNA))\n\nlistap=list(x[\'bcr_patient_barcode\']) #lista pazienti con dati clinici\n\n\nlistanomi=[] #nomi con cui ricavare i dati clinici\nlistanomi01=[ogg_analisi[1]] #nomi con cui ricavare dati di espressione\nfor el in listamiRNA:\n    if el[:-4] in listap:\n        listanomi.append(el[:-4])\n        listanomi01.append(el)\nprint(len(listanomi01))\n        \n#aprire il dataframe di interesse\nd=open_dataframe_gene(ogg_analisi[0],listanomi01, ogg_analisi[2])\n\n\ncolonna1=list(d.loc[gene,])\ncolonna2=list(x1.loc[listanomi, "acronym"])\ncolonna3=list(x1.loc[listanomi, feature])\nprint(len(colonna1),len(colonna2),len(colonna3))\n\n\n\ndf=pd.DataFrame({gene:colonna1, \'tumor\':colonna2, feature:colonna3})\ndf=pulizia_df(df)\nprint(df.shape)\nbox_plot(df)')


# In[ ]:




