#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr as prs
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines.plotting import add_at_risk_counts
from statsmodels.stats import multitest as multi
import numpy as np
import os
import sys
import os
import configparser
from pathlib import Path

# Carica il file di configurazione
config = configparser.ConfigParser()
script_dir = Path(__file__).parent

# Costruisci il percorso relativo al file di configurazione
config_file = script_dir.parent.parent / 'webserver' / 'webserver' / 'conf.ini'

config.read(config_file) 
#directory base
output_data = config['Paths']['output_data']
base_dir = config.get('Paths', 'base_dir', fallback='')

# Costruisci i percorsi completi
def get_full_path(relative_path):
    return os.path.join(base_dir, relative_path)

# In[2]:


def open_dataframe_gene(gene,tumor):
    #aprire il df per i miRNA (isoforme)
    df=pd.read_csv("/mnt/data/notturno/isoforme/mature_Homosapiens.csv",usecols=['miRNA_id'])

    if gene in list(df['miRNA_id']):
        df=pd.read_csv('/mnt/data/notturno/isoforme/Dataframe/df_isoforme_'+tumor+'.csv')
        df=df.set_index('miRNA_region')
        return(df)
    
    #aprire il df per mRNA(gene_expression)
    if gene in open("/mnt/data/notturno/gene_expression/ENSG.txt").read().split("\n"):
        df=pd.read_csv("/mnt/data/notturno/Dataframe_tumorgene/Dataframe_FPKM/Dataframe_FPKM_"+tumor+".csv")
        df=df.set_index("gene_id")
        return (df)
    else: 
        print("per il nome inserito non è disponibile la ricerca")
        return 0


# In[3]:


def converti_name_gene(gene):
    df=pd.read_csv("/mnt/data/notturno/overall_survival_rnainter/gene_symbol_filtrato.csv")
    if gene[:4]!='ENSG':
        df=df.set_index("gene_symbol")
        name=(df.loc[gene])
        return name['ENSG']
    else:
        df=df.set_index("ENSG")
        name=(df.loc[gene])
        return name['gene_symbol']

    
    
def controllo_interazione(mirna, gene):
    #apro df interazioni strong mirna-mrna di rnainter
    df=pd.read_csv("/mnt/data/notturno/overall_survival_rnainter/rnainter_homosapiens_strong.csv")
    
    if mirna in list(df['Interactor1.Symbol']):
        df= df.set_index('Interactor1.Symbol')
        df1=(df.loc[mirna])
        
        for el in list(df1['Interactor2.Symbol']):
            if el==gene_symbol:
                return True
        return False
    else: 
        print('questo mirna non è presente nel file di interazione')
      
    #devo aprire un file he contiene il gene e il mirna che sono presenti nel df


# In[4]:


def dataframe_OStime(dfclinic,tumor):
    df1_mask=dfclinic['type']==tumor
    OS=dfclinic[df1_mask]
    OS=(OS[['bcr_patient_barcode','OS.time']]) 
    OS=OS.set_index('bcr_patient_barcode')
    OS=OS.dropna()
    return(OS)


# In[5]:


def seleziona_samples(dfmrna, dfmirna,listaos,tumor):
    
    #aprire nomi tgca dei mirna (isoforme)
    df=pd.read_csv('/mnt/data/notturno/isoforme/Dataframe/df_isoforme_'+tumor+'.csv')
    listamirna=list(df.columns[1:])
    print("lista campioni miRNA",len(listamirna[1:]))
    
    #nomi dei mrna(gene expression)
    listamrna=open('/mnt/data/notturno/Dataframe_tumorgene/samples/sample_'+tumor+'.txt').read().rstrip().split(",")
    print("lista campioni mrna",len(listamrna[1:]))
    
    #nomi in comune
    lista_corretta=[]
    lista_corretta_01=[]
    for ele in listamrna[1:]:
        if ele[:-4] in listaos:
            if ele in listamirna:
                if ele[-1:]!= "x" and ele[-1:]!="y" and int(ele[-3:-1]) <11:
                    lista_corretta_01.append(ele)
                    lista_corretta.append(ele[:-4])
    return (lista_corretta, lista_corretta_01)
   


# In[6]:


def overall_survival_analysis(mrna, mirna, tumor, dfmrna, dfmirna, OS1, cartella,listaf):
    
    #mrna sopra la media
    i1=dfmrna.loc[mrna,:] > dfmrna.loc[mrna,:].median() #+
    
    #mrna sotto la media
    i2 = dfmrna.loc[mrna,:] < dfmrna.loc[mrna,:].median()  #-
    
    #mirna sopra la media
    j1=dfmirna.loc[mirna,:] > dfmirna.loc[mirna,:].median() #+
    
    #mirna sotto la media
    j2 = dfmirna.loc[mirna,:] < dfmirna.loc[mirna,:].median() #-
    
    
    dfg1=(pd.DataFrame(dict(i1 = i1, i2=i2, j1 = j1, j2 = j2)).reset_index())
    dfg1=dfg1.dropna() 
    
    dfg1['index']=listaf    
    OS=OS1.loc[listaf]
    print(OS.shape)
    
    
    #######################################
    #OS.time del primo gruppo (++)
    d1=(dfg1[(dfg1['i1']==True) & (dfg1['j1']==True)])
    
    OSg1=(OS.loc[list(d1['index'])])
    
    #######################################
    #Os.time del secondo gruppo (+-)
    
    d2=(dfg1[(dfg1['i1']==True) & (dfg1['j2']==True)])
    OSg2=(OS.loc[list(d2['index'])])
    
    
    ##########################################
    #Os.time del terzo gruppo (--)
    
    d3=(dfg1[(dfg1['i2']==True) & (dfg1['j2']==True)])
    OSg3=(OS.loc[list(d3['index'])])
   # display(OSg3['OS.time'])
    
    ##########################################
    #Os.time del quarto gruppo (-+)
    
    d4=(dfg1[(dfg1['i2']==True) & (dfg1['j1']==True)])
    OSg4=(OS.loc[list(d4['index'])])
   # display(OSg4['OS.time'])
    
    
    #grafico 
    k=[OSg1,OSg2, OSg3,OSg4]
    labelmrna=['+','+','-','-']
    labelmirna=['+','-','-','+']
    number=0
    
    f=open(cartella+'result.txt','w')
    for el,OSgroup1 in enumerate(k): 
        
        for ely,OSgroup2 in enumerate(k):
            
            if el !=ely and el<ely:
                
                kmf = KaplanMeierFitter()

                #if np.mean(list(dfmrna.loc[mrna,i2]))>0:
                results = logrank_test((OSgroup1['OS.time']), (OSgroup2['OS.time']), alpha=.95)
                                    
                if results.p_value < 1:
                    number+=1
                    f.write(str(number)+'\tp-value: '+str(results.p_value)+'\n')
                    
                    kmf.fit((OSgroup1), label="mRNA "+labelmrna[el]+" and miRNA "+labelmirna[el])
                    a1 = kmf.plot()

                    kmf.fit((OSgroup2), label="mRNA "+labelmrna[ely]+" and miRNA "+labelmirna[ely])
                    kmf.plot(ax=a1)

                    plt.savefig(cartella+"overallsurvival_"+gene+"_"+tumor+"_"+str(number)+".jpg")
                    plt.clf()
                    
                else:
                    f.write(str(number)+"\tpvalue > 1")
               # else: 
                #    print("MEDIA <0 ??")



#dfclinic=pd.read_csv('/mnt/data/notturno/TCGA-CDR-SupplementalTableS1.csv')
dfclinic=get_full_path(config['clinical']['clinical_OS'])
#input
gene= sys.argv[1]
mirna=sys.argv[2]
tumor=sys.argv[3]

cartella=sys.argv[4]
#cartella='/mnt/data/notturno/web_app/webserver/rolls/static/media/saveanalisi/'+sys.argv[4]+'/' ####da modificare passare come argoemnto dal view.py



if gene[:4]!='ENSG':
    gene_ENSG=converti_name_gene(gene)
    gene_symbol=gene
else:
    gene_symbol=converti_name_gene(gene)
    gene_ENSG=gene
    
#dobbiamo controllare che interagiscano:
if controllo_interazione(mirna,gene_symbol):
    os.mkdir(cartella)
    dfmrna=open_dataframe_gene(gene_ENSG,tumor)
    dfmirna=open_dataframe_gene(mirna,tumor)
    
    OS=dataframe_OStime(dfclinic,tumor)
    listaos=list(OS.index) 
    print("nomi di cui abbiamo os.time",len(listaos)) #nomi senza il 01A finale
    
  #lista pazienti in comune di cui abbiamo i dati di mirna e mrna e OS.time
    listaf, listaf01=(seleziona_samples(dfmrna, dfmirna,listaos,tumor))
    print('pazienti miRNA/mRNA/os: ', len(listaf), len(listaf01))
    
  #OS.time con righe selezionate  
    OS1=OS.loc[listaf,:]
    
  #OVERALL-SURVIVAL
    dfmrna=dfmrna[listaf01]
    overall_survival_analysis(gene_ENSG,mirna,tumor, dfmrna, dfmirna,OS1,cartella,listaf)

else:
    print(mirna,"-",gene,' non interagiscono')
 




