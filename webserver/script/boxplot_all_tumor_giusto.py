#!/usr/bin/env python
# coding: utf-8

# Differential analisys all tumor boxplot

# Passare un gene/miRNA/proteina e una feature

# In[1]:


import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import ranksums
import subprocess
import os
import sys


# In[2]:



def gene_mirna_proteina(gene, cartella):
    
    if gene in open("/mnt/data/notturno/miRNA/namemirna.txt").read().split("\n"):
        os.mkdir(cartella)
        return(['miRNA',"miRNA_ID",'/mnt/data/notturno/miRNA/DataFrameTCGA_miRNA.csv'])
    
    if gene in open("/mnt/data/notturno/gene_expression/ENSG.txt").read().split("\n"):
        os.mkdir(cartella)
        listageni=open("/mnt/data/notturno/gene_expression/ENSG.txt").read().strip().split("\n")
        posizione="-e "+str(listageni.index(gene)+1)+"p"
        
        #creiamo un dataframe piu piccolo dove c'è solo la riga del gene che è stato selezionato
       
        path=cartella+'/'+gene+"_df.txt"
        out_file=open(path,"w")
        subprocess.call(["sed","-n", "-e 1p", posizione, "/mnt/data/notturno/gene_expression/DataFrameTCGA_gene.csv"],stdout=out_file)
        
        return(['gene','gene_id',path])
    
    
    if gene in open("/mnt/data/notturno/protein/namepeptide.csv").read().split("\n"):
        os.mkdir(cartella)
        return(['protein', "peptide_target",'/mnt/data/notturno/protein/DataFrameTCGA_protein.csv'])
    
    else:
        return(0) #la ricerca non è disponibile per il nome inserito


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


# In[4]:


def box_plot(df1, cartella, gene, feature):
    
    sns.set(rc={'figure.figsize':(25.7,8.27)})
    sns.set_style("white")

    my_order = df1.groupby(by=["tumor"])[gene].median().iloc[::-1].sort_values().index

    ax=sns.boxplot(x="tumor", y=gene, hue=feature, data=df1, palette="Set2", width=0.7, order=my_order)
    
    ax.set_yscale("log")
    
    plt.savefig(cartella+'/'+gene+'_'+feature+'.jpg')


# In[5]:

def df_feature_age(x):
    if feature=="age_at_initial_pathologic_diagnosis":
        df=pd.read_csv('/mnt/data/notturno/mediana_age_tumor.csv')
        df=df.set_index('acronym')
        
        lista_age=[]
        x1=x
        x1 = x1.drop(x1[x1[feature]== '[Not Available]'].index)

        for tumor,age in zip(x1['acronym'],x1[feature]):
            age=int(age)
            mediana=df.loc[tumor,'median_age']
            
            if age<=mediana:
                lista_age.append("under")
            else:
                lista_age.append("over")

        x1=x1.rename(columns={"age_at_initial_pathologic_diagnosis":"age"})
        x1[feature]=lista_age

        return(x1) 


#per la feature patient status
def crealista_patient(dffeat,listamiRNA,ogg):
    lista=list(dffeat.index)
        
    listaf01=[ogg[1]]
    listaf=[]
    for ele in listamiRNA:
        if ele[:-4] in lista and ele[-1]!="x" and ele[-1]!='y':
            listaf01.append(ele)
            listaf.append(ele[:-4])
            
    
    return listaf01,listaf

def pulisci_df(df):
    for tumor in set(list(df.tumor)):
        df1=df
        df1_mask=df['tumor'] == tumor
        df1=df1[df1_mask]

        lista_feature=((list(df1[feature])))
        t=0
        c=0
        for ele in lista_feature:
            if ele== 'Tumor':
                t+=1
            else:
                c+=1       

        #print(tumor, t, c)
        if t==0 or c==0:
            df.drop(df.index[df['tumor'] == tumor], inplace = True)
            
    return df

#ranksum test per p-value
def p_value(df, cartella):
    f=open(cartella+'/result.txt','w')
    f.write('tumor'+'\t'+'p-value'+'\n')
    tumori=set(list(df.tumor))
    for tumor in tumori:
        p=list(set(df[feature]))

        df1=df
        df1_mask=df1['tumor'] == tumor
        df1=df1[df1_mask]


        filtered_df = df1[df1[feature].eq(p[0])]
        list0= list(filtered_df[gene])


        filtered_df = df1[df1[feature].eq(p[1])]
        list1= list(filtered_df[gene])

        w, p = ranksums(list0, list1)
        print(tumor, p)
        f.write(tumor+"\t"+str(p)+"\n")
        


# In[8]:

path="/mnt/data/notturno/"

gene= sys.argv[1]
feature= sys.argv[2]

cartella='/mnt/data/notturno/web_app/webserver/rolls/static/media/saveanalisi/'+sys.argv[3]


x=pd.read_csv("/mnt/data/notturno/clinical_PANCAN_patient_with_followup_modificato.csv")
x=x.set_index("bcr_patient_barcode")


if feature=="age_at_initial_pathologic_diagnosis":
    x1=df_feature_age(x)
else:
    x1=x


#aprire df
ogg_analisi=gene_mirna_proteina(gene, cartella) #es. miRNA, miRNA_ID, 'pathdataframe.csv'

if ogg_analisi!=0:

    #lista samples df expression gene/miRNA/prot
    dd=pd.read_csv(ogg_analisi[2],nrows=1)
    listamiRNA=list(dd.columns[1:])


    #patient status
    if feature== 'patient_status':    
        listanomi01, listanomi= crealista_patient(x, listamiRNA,ogg_analisi)
        
        #aprire df espressione
        d=open_dataframe_gene(ogg_analisi[0],listanomi01, ogg_analisi[2])
        

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
        d=open_dataframe_gene(ogg_analisi[0],listanomi01, ogg_analisi[2])



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
        df=pulisci_df(df)  



    box_plot(df, cartella, gene, feature)


    p_value(df, cartella)

else:
    print('Ricerca non disponibile per il nome inserito')