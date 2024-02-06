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

# In[1]:


import pandas as pd
import os 
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ranksums
import sys


# In[2]:


def open_dataframe_gene(gene,tumor):
    #controllo df geni
    if gene[:4]=='ENSG':
        listageni=open("/mnt/data/notturno/gene_expression/ENSG.txt").read().split("\n")
        c=0
        for g in listageni:
            if gene in g:
                gene=g
                c=1
        if c==1:
            if feature == 'patient_status':
                posizione="-e "+str(listageni.index(gene)+1)+"p"
                #creiamo un dataframe piu piccolo dove c'è solo la riga del gene che è stato selezionato

                path=cartella+'/'+gene+"_df.txt"
                out_file=open(path,"w")
                subprocess.call(["sed","-n", "-e 1p", posizione, "/mnt/data/notturno/gene_expression/DataFrameTCGA_gene.csv"],stdout=out_file)

                df=pd.read_csv(path)
                df=df.set_index("gene_id")
                return(df, 'gene', gene)
            else:
                df=pd.read_csv("/mnt/data/notturno/Dataframe_tumorgene/Dataframe_FPKM/Dataframe_FPKM_"+tumor+".csv")
                df=df.set_index("gene_id")
                return (df, 'gene', gene) 

    


    #controllo aprire df miRNA
    if gene in open("/mnt/data/notturno/miRNA/namemirna.txt").read().split("\n"):
            df=pd.read_csv('/mnt/data/notturno/miRNA/DataFrameTCGA_miRNA.csv')
            df=df.set_index('miRNA_ID')
            return(df, 'miRNA', gene)

    #controllo df proteine
    if gene in open("/mnt/data/notturno/protein/namepeptide.csv").read().split("\n"):
        df=pd.read_csv("/mnt/data/notturno/protein/DataFrameTCGA_protein.csv")
        df=df.set_index('peptide_target')
        return (df, 'protein', gene)
    else: 
        #print("per il nome inserito non è disponibile la ricerca")
        return 0


# In[3]:


def df_feature(ogg, df, tumor, feature):
   
    if ogg =='miRNA' or ogg=='protein' or ogg=='gene':
        x=pd.read_csv("/mnt/data/notturno/clinical_PANCAN_patient_with_followup_modificato.csv")
        df1_mask=x['acronym']== tumor
        dfclinical=x[df1_mask]

        if feature=="age_at_initial_pathologic_diagnosis":
            median=dfclinical[feature].median()
            
            lista_age=[]
            for ele in dfclinical[feature]:
                ele=int(ele)
                if ele<=median:
                    lista_age.append("under")
                else:
                    lista_age.append("over")
                    
            dfclinical=dfclinical.rename(columns={"age_at_initial_pathologic_diagnosis":"age"})
            dfclinical[feature]=lista_age

        return (dfclinical)
        


# In[4]:


def crealista(dffeat,df, feature):
    lista=list(dffeat['bcr_patient_barcode'])
    
    listaf01=[]
    listaf=[]
    
    for ele in df.columns:
        if ele[:-4] in lista and ele[-1]!="x" and ele[-1]!='y':
            if int(ele[-3:-1])<11:
                if ele[:-4] not in listaf:
                    listaf.append(ele[:-4])
                    listaf01.append(ele)

    if len(listaf)==len(listaf01):
        return listaf, listaf01
    else:
        return 0

#per la feature patient status
def crealista_patient(dffeat,df):
    lista=list(dffeat['bcr_patient_barcode'])
        
    listaf01=[]
    listaf=[]
    for ele in df.columns:
        if ele[:-4] in lista and ele[-1]!="x" and ele[-1]!='y':
            listaf01.append(ele)

    
    return listaf01
   
#funzioni per interattori dei miRNA#######################

def lista_interattori(mirna, df):
    #apro df interazioni strong mirna-mrna di rnainter
    
    if mirna in list(df['Interactor1.Symbol']):
        df1=(df[df['Interactor1.Symbol']==mirna])
        #display(df1)
        interactor=(list(df1['Interactor2.Symbol']))
        if (len(interactor))>1:
            return interactor
        else:
            return 0
    else:
        return 0
    
def conversione_ensg(lista_gs):
    df=pd.read_csv("/mnt/data/notturno/overall_survival_rnainter/gene_symbol_filtrato.csv")
    lista_E=(open("/mnt/data/notturno/gene_expression/ENSG.txt").read().split('\n'))
    name_ENSG=[]
    ld=[]
    for gene in df['gene_symbol']:
        if gene in lista_gs:
            if gene not in ld:
                ld.append(gene)
                #il gene ensg deve essere presente nella lista di geni con dati di espressione
                df1=(df[df['gene_symbol']==gene])
               
                g=list(df1['ENSG'])

                for i in g:
                    if i in lista_E:
                        name_ENSG.append(i)
                
               
     
    
    return name_ENSG
##########################################

# In[5]:


#parametri da inserire:
gene= sys.argv[1]
tumor= sys.argv[2]
feature= sys.argv[3]

cartella='/mnt/data/notturno/web_app/webserver/rolls/static/media/saveanalisi/'+ sys.argv[4]
os.mkdir(cartella)

#df di interesse:
ogg=open_dataframe_gene(gene,tumor)
gene=ogg[2]

if ogg!=0:
    df=ogg[0]
    dffeat=df_feature(ogg[1], df, tumor, feature) #df features
    
    if feature== 'patient_status':
        listaf0= crealista_patient(dffeat, df)
        
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
        listaf, listaf01= crealista(dffeat, df, feature)#lista di samples di cui abbiamo le features e dati di espressione
        dffeat=dffeat[['bcr_patient_barcode',feature]]

        #creo un dframe finale
        dframe=pd.DataFrame({'nome01':listaf01, 'bcr_patient_barcode':listaf})
        d= dframe.merge(dffeat)

        lista_exp=[] #lista valori espressione
        for ele in d['nome01']:
            lista_exp.append(df.loc[gene,ele])

        d[gene]=lista_exp

    #plt
    ax=sns.boxplot(x=feature, y=gene, data=d, width=0.7,hue=feature)
    if ogg[1]== "miRNA":
        ax.set_yscale("log")

    plt.savefig(cartella+'/'+gene+'_'+feature+'.jpg')


    #ranksum test per p-value
    p=list(set(d[feature]))
    
    df1_mask=d[feature]== p[0]
    dp0=d[df1_mask]
    
    df1_mask=d[feature]== p[1]
    dp1=d[df1_mask]

    w, p = ranksums(list(dp0[gene]), list(dp1[gene]))

    print(p)

else: 
    print(0)

if ogg[1]=='miRNA':
        #df con interazioni mirna(actor1)-mrna(actor2)
        df_interaction=pd.read_csv('/mnt/data/notturno/overall_survival_rnainter/rnainter_homosapiens_strong.csv')

        #dato un mirna cerchiamo tutte le sue isoforme
        isoforme=[]
        mirna=gene
        
        lista=list(set(list(df_interaction['Interactor1.Symbol'])))
        for ele in lista[1:]:
            if mirna in ele:
                isoforme.append(ele)

        if len(isoforme)!=0:
            #cerchiamo interattori del miRNA

            #dataframe espressione FPKM
            df_exp=pd.read_csv('/mnt/data/notturno/Dataframe_tumorgene/Dataframe_FPKM/Dataframe_FPKM_'+tumor+'.csv')
            df_exp=df_exp.set_index('gene_id')
            
            for mirna in isoforme:
                #print(mirna)
                lista_interact= lista_interattori(mirna, df_interaction)
                #print(lista_interact)
                if lista_interact!=0:
                    lista_ENSG=conversione_ensg(lista_interact)


                    df_finale=df_exp.loc[lista_ENSG,:]

                    #rimuovo i geni che hanno espressione per tutti i campioni =0
                    df0=df_finale.loc[(df_finale != 0).any(axis=1)]

                    df0.to_csv(cartella+'/result.txt')

                    df=df0.T
                    df=df.dropna()
                    #df.to_csv(cartella+'/'+mirna+'_dfexpression.csv')
                    sns.clustermap(df, standard_scale=1, method="ward", cmap="mako",row_cluster=False, robust=True)
                    plt.savefig(cartella+'/'+mirna+'_cluestermap.jpg')
                   
               