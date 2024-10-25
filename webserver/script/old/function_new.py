import pandas as pd
import os 
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import subprocess
from pathlib import Path
import configparser
from scipy.stats import ranksums
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines.plotting import add_at_risk_counts
from statsmodels.stats import multitest as multi

#######################################################################################
#                              Carica il file di configurazione                       #
#######################################################################################

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

# Accesso ai valori e costruzione dei percorsi
if 'miRNA' in config:
    miRNA_name = get_full_path(config['miRNA']['name'])
    miRNA_dataframe = get_full_path(config['miRNA']['dataframe'])

if 'gene' in config:
    gene_name_ENSG = get_full_path(config['gene']['name_ENSG'])
    gene_dataframe = get_full_path(config['gene']['dataframe'])
    gene_dataframe_FPKM = get_full_path(config['gene']['dataframe_FPKM'])

if 'protein' in config:
    protein_name = get_full_path(config['protein']['name'])
    protein_dataframe = get_full_path(config['protein']['dataframe'])

if 'clinical' in config:
    clinical_data= get_full_path(config['clinical']['dati_clinici'])



#######################################################################################
#                                      FUNCTION                                       #
#######################################################################################

#######  ->                           Boxplot_all_tumor                     <-  #######
 
def read_clinical_data():
    x=pd.read_csv(clinical_data)
    x=x.set_index("bcr_patient_barcode")
    return(x)

def detect_if_gene_mirna_proteina(gene, cartella):
    
    if gene in open(miRNA_name).read().split("\n"):
        os.mkdir(cartella)
        return(['miRNA',"miRNA_ID",miRNA_dataframe])
    
    if gene in open(gene_name_ENSG).read().split("\n"):
        os.mkdir(cartella)
        listageni=open(gene_name_ENSG).read().strip().split("\n")
        posizione="-e "+str(listageni.index(gene)+1)+"p"
        
        #creiamo un dataframe piu piccolo dove c'è solo la riga del gene che è stato selezionato
        path=cartella+'/'+gene+"_df.txt"
        out_file=open(path,"w")
        subprocess.call(["sed","-n", "-e 1p", posizione, gene_dataframe ],stdout=out_file)
        
        return(['gene','gene_id',path])

    if gene in open(protein_name).read().split("\n"):
        os.mkdir(cartella)
        return(['protein', "peptide_target",protein_dataframe])
    
    else:
        return(0) #la ricerca non è disponibile per il nome inserito

def open_dataframe_gene_boxplot_all_tumor(gene,listanomi01, path_dataframe,index):
    if gene == 'miRNA':
        df=pd.read_csv(path_dataframe,usecols=listanomi01)
        df=df.set_index(index)
        return(df)
    if gene == 'gene':
        df=pd.read_csv(path_dataframe,usecols=listanomi01)
        df=df.set_index(index)
        return (df)
    if gene =='protein':
        df=pd.read_csv(path_dataframe,usecols=listanomi01)
        df=df.set_index(index)
        return (df)
    else: 
        print("per il nome inserito non è disponibile la ricerca")
        return 0



def box_plot_all_tumor(df1, cartella, gene, feature):
    
    sns.set(rc={'figure.figsize':(25.7,8.27)})
    sns.set_style("white")

    my_order = df1.groupby(by=["tumor"])[gene].median().iloc[::-1].sort_values().index
    ax=sns.boxplot(x="tumor", y=gene, hue=feature, data=df1, palette="Set2", width=0.7, order=my_order)
    ax.set_yscale("log")
    
    plt.savefig(cartella+'/'+gene+'_'+feature+'.jpg')


def df_feature_age(x,feature):
    if feature=="age_at_initial_pathologic_diagnosis":
        df=pd.read_csv('mediana_age_tumor.csv')
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
def crealista_patient_status(dffeat,listamiRNA,ogg):
    lista=list(dffeat.index)
        
    listaf01=[ogg[1]]
    listaf=[]
    for ele in listamiRNA:
        if ele[:-4] in lista and ele[-1]!="x" and ele[-1]!='y':
            listaf01.append(ele)
            listaf.append(ele[:-4])
            
    return listaf01,listaf


def pulisci_df(df,feature):
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
def p_value(df, cartella,feature,gene):
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


#######################################################################################

#######  ->                       Boxplot_single_tumor                      <-  #######
def what_is_my_object_gene(gene):
    if 'ENSG' in gene:
        print('ENSG',gene)
        df_ensg= pd.read_csv(gene_name_ENSG,sep='\t')
        print(df_ensg)

        result_index = df_ensg[(df_ensg['gene_id_version'] == gene) | (df_ensg['gene_id'] == gene)].index
        if not result_index.empty:
            gene_version=(df_ensg.loc[result_index[0],'gene_id_version'])
            return (gene_version,'gene','gene_id')

       
    if gene in open(protein_name).read().split("\n"):
        return(gene,'protein','peptide_target')

    if gene in open(miRNA_name).read().split("\n"):
        return (gene,'miRNA','miRNA_ID')

    else:
        return(0)

def open_dataframe(gene,tumor,feature,cartella):
    input=what_is_my_object_gene(gene)
    if input!=0:
        if input[1]=='gene':
            if feature == 'patient_status':
                listageni=open(gene_name_ENSG).read().strip().split("\n")
                posizione="-e "+str(listageni.index(gene)+1)+"p"

                #creiamo un dataframe piu piccolo dove c'è solo la riga del gene che è stato selezionato

                path=cartella+'/'+gene+"_df.txt"
                out_file=open(path,"w")
                subprocess.call(["sed","-n", "-e 1p", posizione, gene_dataframe],stdout=out_file)

                df=pd.read_csv(path)
                df=df.set_index("gene_id")
                return(df, 'gene')
            else:
                path_df="Dataframe_FPKM_"+tumor+".csv"
                df=pd.read_csv(os.path.join(gene_dataframe_FPKM,path_df))
                df=df.set_index("gene_id")
                return (df, 'gene') 
        else:
            df=df=pd.read_csv(input[3])
            df=df.set_index(input[2])
            return(df,input[1])
        
    else: 
        print("non è disponibile la ricerca per il nome inserito")
        return(0)
    


def open_dataframe_gene(gene,tumor,feature,cartella):

    if gene in open(miRNA_name).read().split("\n"):
        df=pd.read_csv(miRNA_dataframe)
        df=df.set_index('miRNA_ID')
        return(df, 'miRNA')

    if gene in open(gene_name_ENSG).read().split("\n"):
        if feature == 'patient_status':
            listageni=open(gene_name_ENSG).read().strip().split("\n")
            posizione="-e "+str(listageni.index(gene)+1)+"p"

            #creiamo un dataframe piu piccolo dove c'è solo la riga del gene che è stato selezionato

            path=cartella+'/'+gene+"_df.txt"
            out_file=open(path,"w")
            subprocess.call(["sed","-n", "-e 1p", posizione, gene_dataframe],stdout=out_file)

            df=pd.read_csv(path)
            df=df.set_index("gene_id")
            return(df, 'gene')
        else:
            path_df="Dataframe_FPKM_"+tumor+".csv"
            df=pd.read_csv(os.path.join(gene_dataframe_FPKM,path_df))
            df=df.set_index("gene_id")
            return (df, 'gene') 

    if gene in open(protein_name).read().split("\n"):
        df=pd.read_csv(protein_dataframe)
        df=df.set_index('peptide_target')
        return (df, 'protein')
    else: 
        print("per il nome inserito non è disponibile la ricerca")
        return 0
    


def df_feature(ogg, tumor, feature):
   
    if ogg =='miRNA' or ogg=='protein' or ogg=='gene':
        x=pd.read_csv(clinical_data)
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

#######################################################################################

#######  ->                        Overall Survival                         <-  #######

def open_dataframe_gene_overall(gene,tumor):
    if gene in open(miRNA_name).read().split("\n"):
        df=pd.read_csv(miRNA_dataframe)
        df=df.set_index('miRNA_ID')
        return(df)
    if gene in open(gene_name_ENSG).read().split("\n"):
        gene_dataframe_FPKM_tumor=gene_dataframe_FPKM+"_"+tumor+".csv"
        df=pd.read_csv(gene_dataframe_FPKM_tumor)
        #df=pd.read_csv("/mnt/data/notturno/Dataframe_tumorgene/Dataframe_FPKM/Dataframe_FPKM_"+tumor+".csv")
        df=df.set_index("gene_id")
        return (df)
    if gene in open(protein_name).read().split("\n"):
        df=pd.read_csv(protein_dataframe)
        df=df.set_index('peptide_target')
        return (df)
    else: 
        print("per il nome inserito non è disponibile la ricerca")
        return 0



def dataframe_OStime(tumor):
    dfclinic=pd.read_csv(clinical_data)    
    OS=(dfclinic[['bcr_patient_barcode','OS.time']]) 
    df1_mask=dfclinic['type']==tumor
    OS=dfclinic[df1_mask]
    OS=(OS[['bcr_patient_barcode','OS.time']]) 
    OS=OS.set_index('bcr_patient_barcode')
    #display(OS)
    OS=OS.dropna()
    return(OS)


def overall_survival_analysis(m,tumor,feature,cartella,df1,OS1,gene):
    
    
    i1=df1.loc[m,:] > df1.loc[m,:].median()
    i2 = df1.loc[m,:] < df1.loc[m,:].median() 
    
    kmf = KaplanMeierFitter()
    

    #if np.mean(list(df1.loc[m,i2]))>0:
    results = logrank_test((OS1[i1]), (OS1[i2]),list(df1.loc[m,i1]),list(df1.loc[m,i2]), alpha=.95)
    
    if results.p_value < 1:
        os.mkdir("../rolls/static/media/saveanalisi/"+cartella)
        print("p-value:",results.p_value)
        
        kmf.fit((OS1[i1]), list(df1.loc[m,i1]), label="Higher expression")
        a1 = kmf.plot()

        kmf.fit((OS1[i2]),list(df1.loc[m,i2]) , label="Lower expression")
        kmf.plot(ax=a1)
        plt.savefig("../rolls/static/media/saveanalisi/"+cartella+"/overallsurvival_"+gene+"_"+tumor+".png")
        
    else:
        print("ANALISI NON VALIDA pvalue>1")
  

  #######################################################################################