import pandas as pd
import os 
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import subprocess


def open_dataframe_gene(gene,tumor,feature,cartella):

    if gene in open("miRNA/namemirna.txt").read().split("\n"):
        df=pd.read_csv('miRNA/DataFrameTCGA_miRNA.csv')
        df=df.set_index('miRNA_ID')
        return(df, 'miRNA')

    if gene in open("gene_expression/ENSG.txt").read().split("\n"):
        if feature == 'patient_status':
            listageni=open("gene_expression/ENSG.txt").read().strip().split("\n")
            posizione="-e "+str(listageni.index(gene)+1)+"p"

            #creiamo un dataframe piu piccolo dove c'è solo la riga del gene che è stato selezionato

            path=cartella+'/'+gene+"_df.txt"
            out_file=open(path,"w")
            subprocess.call(["sed","-n", "-e 1p", posizione, "gene_expression/DataFrameTCGA_gene.csv"],stdout=out_file)

            df=pd.read_csv(path)
            df=df.set_index("gene_id")
            return(df, 'gene')
        else:
            df=pd.read_csv("Dataframe_tumorgene/Dataframe_FPKM/Dataframe_FPKM_"+tumor+".csv")
            df=df.set_index("gene_id")
            return (df, 'gene') 

    if gene in open("protein/namepeptide.csv").read().split("\n"):
        df=pd.read_csv("protein/DataFrameTCGA_protein.csv")
        df=df.set_index('peptide_target')
        return (df, 'protein')
    else: 
        print("per il nome inserito non è disponibile la ricerca")
        return 0
    


def df_feature(ogg, tumor, feature):
   
    if ogg =='miRNA' or ogg=='protein' or ogg=='gene':
        x=pd.read_csv("clinical_PANCAN_patient_with_followup_modificato.csv")
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
    

def crealista(dffeat,df,feature):
    #lista di samples di cui abbiamo le features e dati di espressione
    
    #caso specifico se abbiamo feature:patient_status
    if feature=='patient_status':
        lista=list(dffeat['bcr_patient_barcode'])
        listaf01=[]
        listaf=[]
        for ele in df.columns:
            if ele[:-4] in lista and ele[-1]!="x" and ele[-1]!='y':
                listaf01.append(ele)

        
        return listaf01

    #altrimenti con tutte le altre feature:
    else:
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


