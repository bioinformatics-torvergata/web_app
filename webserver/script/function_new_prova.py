import pandas as pd
import os 
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import subprocess
from pathlib import Path
import configparser
from scipy.stats import ranksums
import numpy as np
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines.plotting import add_at_risk_counts
from statsmodels.stats import multitest as multi
import plotly.express as px
import plotly.graph_objects as go 
import shutil
#from statannotations.Annotator import Annotator
from decimal import Decimal

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
    mapped_ENSG=get_full_path(config['gene']['mapped_ENSG'])

if 'protein' in config:
    protein_name = get_full_path(config['protein']['name'])
    protein_dataframe = get_full_path(config['protein']['dataframe'])

if 'clinical' in config:
    clinical_data= get_full_path(config['clinical']['dati_clinici'])
    dati_age= get_full_path(config['clinical']['dati_age'])
    clinical_OS= get_full_path(config['clinical']['clinical_OS'])

if 'os' in config:
    os_pathway=get_full_path(config['os']['os_pathway'])

if 'TCGA' in config:
    TCGA_path=get_full_path(config['tcga']['split_count'])



#######################################################################################
#                                      FUNCTION                                       #
#######################################################################################

#######  ->                           Boxplot_all_tumor                     <-  #######
 
def read_clinical_data():
    x=pd.read_csv(clinical_data)
    x=x.set_index("bcr_patient_barcode")
    return(x)



def gene_mirna_proteina(gene, cartella,control):
    if control == 'protein':
        with open(protein_name) as f:
            protein_names = f.read().splitlines()
        for protein in protein_names:
            if protein.lower() == gene.lower():  # Use case-insensitive comparison to find if gene is in the list
                os.mkdir(cartella)
                return ['protein', "peptide_target", protein_dataframe, protein] # Return the correct-cased protein name along with other info
        return 0   # If gene not found in any case, return 0
    
    else:
        #controllare se è gene symbol:
        df_ensg= pd.read_csv(gene_name_ENSG,sep='\t')

        result_index = df_ensg[(df_ensg['gene_id_version'].str.lower() == gene.lower()) | (df_ensg['gene_id'].str.lower() == gene.lower()) | (df_ensg['gene_symbol'].str.lower() == gene.lower())].index
    
        #result_index = df_ensg[(df_ensg['gene_id_version'] == gene) | (df_ensg['gene_id'] == gene) | (df_ensg['gene_symbol'] == gene)].index
        if not result_index.empty:
            gene_version=df_ensg.loc[result_index[0],'gene_id_version']
            indice=int(result_index[0])
            os.mkdir(cartella)
            posizione="-e "+str(indice+2)+"p"
                #creiamo un dataframe piu piccolo dove c'è solo la riga del gene che è stato selezionato
                
            path=cartella+'/'+gene+"_df.txt"
            out_file=open(path,"w")
            subprocess.call(["sed","-n", "-e 1p", posizione,gene_dataframe],stdout=out_file)
            return(['gene','gene_id',path,gene_version])   


        #controllo df miRNA
        with open(miRNA_name) as f:
            miRNA_names = f.read().splitlines()
            for mirna in miRNA_names:
                if mirna.lower() == gene.lower():
                    os.mkdir(cartella)
                    return(['miRNA',"miRNA_ID",miRNA_dataframe,mirna])
            else:
                return(0) #la ricerca non è disponibile per il nome inserito




def open_dataframe_gene_boxplot_all_tumor(gene,listanomi01, path_dataframe, index):
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
        #print("per il nome inserito non è disponibile la ricerca")
        return 0



def box_plot_all_tumor(df1, cartella, gene, feature,type_gene,gene_input):
    
    sns.set_theme(rc={'figure.figsize':(25.7,8.27)})
    sns.set_style("white")

    my_order = df1.groupby(by=["tumor"])[gene].median().iloc[::-1].sort_values().index
    ax=sns.boxplot(x="tumor", y=gene, hue=feature, data=df1, palette="Set2", width=0.7, order=my_order)
    if type_gene == 'miRNA':
        ax.set_yscale("log")
    
    ax.set_ylabel(f"{gene_input}", fontsize=14) 
    plt.savefig(cartella+'/'+gene_input+'_'+feature+'.jpg',dpi=300)


def df_feature_age(x,feature):
    if feature=="age_at_initial_pathologic_diagnosis":
        df=pd.read_csv(dati_age)
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

        w, p_value = ranksums(list0, list1)
        print(tumor, p)
        #f.write(tumor+"\t"+str(p_value)+"\n")
        p_value_scientific = "%.2E" %Decimal(p_value)
        
        # print(tumor, p_value_scientific)
        f.write(tumor + "\t" + p_value_scientific + "\n")
       


#######################################################################################

#######  ->                       Boxplot_single_tumor                      <-  #######

def what_is_my_object_gene(gene,control):
    if control=='protein':
        with open(protein_name) as f:
            protein_names = f.read().splitlines()
        for protein in protein_names:
            if protein.lower() == gene.lower():
                return(protein,'protein','peptide_target',protein_dataframe)
        else:
            return(0) #nome errato
   
    else:
         #implementato per prendere in input anche l'ENSG inserito senza versione.
        df_ensg= pd.read_csv(gene_name_ENSG,sep='\t')
        result_index = df_ensg[(df_ensg['gene_id_version'].str.lower() == gene.lower()) | (df_ensg['gene_id'].str.lower() == gene.lower()) | (df_ensg['gene_symbol'].str.lower() == gene.lower())].index
        #result_index = df_ensg[(df_ensg['gene_id_version'] == gene) | (df_ensg['gene_id'] == gene) | (df_ensg['gene_symbol'] == gene)].index
        if not result_index.empty:
            gene_version=df_ensg.loc[result_index[0],'gene_id_version']
            return (gene_version,'gene','gene_id',int(result_index[0]))
    
       #controllo df miRNA
        with open(miRNA_name) as f:
            miRNA_names = f.read().splitlines()
            for mirna in miRNA_names:
                if mirna.lower() == gene.lower():
                    return (mirna,'miRNA','miRNA_ID',miRNA_dataframe)
            else:
                return(0) #nome errato
    
def open_df_gene(input,tumor,feature,cartella):
    if feature == 'patient_status':
        #print(input)
        posizione="-e "+str(input[3]+2)+"p"
        #creiamo un dataframe piu piccolo dove c'è solo la riga del gene che è stato selezionato
        path=cartella+'/'+str(input[0])+"_df.txt"
        out_file=open(path,"w")
        subprocess.call(["sed","-n", "-e 1p", posizione, gene_dataframe],stdout=out_file)
        return(path)
                
    else:
        path_df="Dataframe_FPKM_"+tumor+".csv"
        path=os.path.join(gene_dataframe_FPKM,path_df)
        return(path)


def open_dataframe(gene,tumor,feature,cartella,control):
    input=what_is_my_object_gene(gene,control)
    if input!=0:
        if input[1]=='gene':
            #print(feature)
            path=open_df_gene(input,tumor,feature,cartella) #bisogna vedere in base al tipo di feature che viene passata
            df=pd.read_csv(path)
            df=df.set_index(input[2])
            return(df, input[1],input[0])
        else:
            df=df=pd.read_csv(input[3])
            df=df.set_index(input[2])
            return(df,input[1],input[0])
    else: 
        # print("not avaible name for analysis")
        return(0)
    

def plotly_plot(feature,d, gene,cartella,ogg,gene_input):
        fig = px.scatter(x=range(10), y=range(10))
        fig=px.box(d,y=gene,x=feature,color=feature) #points = 'all'
        if ogg[1]== "miRNA":
                fig.update_layout(yaxis_type="log")
        fig.update_layout(yaxis_title=gene_input)    
        fig.write_html(cartella+'/'+gene_input+'_'+feature+'.html')



def ranksum_test(gene,d,feature,cartella,tumor):
    #ranksum test per p-value
    p=list(set(d[feature]))
    
    df1_mask=d[feature]== p[0]
    dp0=d[df1_mask]
    dp0[gene]=dp0[gene].fillna(0)


    df1_mask=d[feature]== p[1]
    dp1=d[df1_mask]
    dp1[gene]=dp1[gene].fillna(0)

    w, p_value = ranksums(list(dp0[gene]), list(dp1[gene]))
    if p_value=='':
        p_value='not significant'
    
    # Dimensione dei campioni
    size_dp0 = len(dp0)
    values_dp0=(list(dp0[gene]))
    size_dp1 = len(dp1)
    values_dp1=(list(dp1[gene]))

    if all(value == 0.0 for value in values_dp0) and all(value == 0.0 for value in values_dp1):
        return(2)

    else:
        p_value_scientific = '%.2E' %Decimal(p_value)
        # Creazione del dizionario dei risultati
        results = {
            "Gene": gene,
            "Tumor":tumor,
            "Feature": feature,
            "Group 1":p[0],
            "Group 2":p[1],
            "n.Sample Group 1": size_dp0,
            "n.Sample Group 2": size_dp1,
            "P-value": p_value_scientific
        }
        
        # Creazione di un DataFrame per facilitare l'output
        results_df = pd.DataFrame([results])
        
        # Scrittura su file CSV
        output_file=os.path.join(cartella,'Result.txt')
        results_df.to_csv(output_file, index=False)
        
        return(p_value)





def df_feature(ogg, tumor, feature):
    x=pd.read_csv(clinical_data)
    df1_mask=x['acronym']== tumor
    dfclinical=x[df1_mask]

    if feature=="age_at_initial_pathologic_diagnosis":
       # dfclinical = dfclinical.dropna(subset=['age_at_initial_pathologic_diagnosis'])
        dfclinical[feature] = pd.to_numeric(dfclinical[feature], errors='coerce').astype('Int64')


        median=dfclinical[feature].median()
        
        lista_age=[]
        for ele in dfclinical[feature]:
            
            ele=int(ele)
            if ele<=median:
                lista_age.append("under")
            else:
                lista_age.append("over")
                
        dfclinical=dfclinical.rename(columns={feature:"age"})
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

#######  ->                               Deseq2                            <-  #######

def copyfile(tumor,pathfiles,dir_saveresults):
    # print('siamo nel copy file')
    # print(' i file li pesca da qui ', pathfiles)
    files=os.listdir(pathfiles)
    print(files)
    for file in files:
        
        if tumor in file:
            
            path_file=os.path.join(pathfiles,file)
            # print('prende da qui : ',path_file)
            
            copy_filepath=os.path.join(dir_saveresults, os.path.basename(file))
            # print('copia qui: ',copy_filepath)
       
            shutil.copy(path_file,copy_filepath )

    #copy file zip:
    result_zip=os.path.join(pathfiles,'results.zip')
    subprocess.run(["cp", result_zip, dir_saveresults])#,cwd=dir) 

def plotly_volcano(df,cartella,tumor):

    significance_threshold = -np.log10(0.05)
    fold_change_threshold = 1
    df['color'] = 'grey'  # default colore
    df.loc[(df['padj'] > significance_threshold), 'color'] = 'grey'  
    df.loc[(df['padj'] < significance_threshold), 'color'] = 'grey'  
    df.loc[(df['log2FoldChange'] >= fold_change_threshold) & (df['padj'] > significance_threshold), 'color'] = 'red' # Up-regolati
    df.loc[(df['log2FoldChange'] <= -fold_change_threshold) & (df['padj'] > significance_threshold), 'color'] = 'blue'  # Down-regolati
   
    fig=go.Figure()
    trace1=go.Scatter(
        x=df['log2FoldChange'],
        y=df['padj'],
        mode='markers',
        hovertext=list(df['GeneSymbol']),
        marker=dict(
            color=df['color'],  # Usa la colonna color per colorare i punti
            size=10
        ),
    )
    fig.add_trace(trace1)
    fig.add_hline(y=(significance_threshold),line_dash="dash")# Padj cutoff
    fig.add_vline(x=fold_change_threshold, line_dash="dash")# line_color="black")  # Fold change positivo cutoff
    fig.add_vline(x=-fold_change_threshold, line_dash="dash")#, line_color="black")  # Fold change negativo cutoff

    fig.update_layout(
        title="Volcano Plot",
        xaxis_title="log2 Fold Change",
        yaxis_title="-log10(padj)"
    )
    fig.write_html(cartella+'/'+tumor+'.html')



    


#######################################################################################

#######  ->                        Overall Survival                         <-  #######

def open_dataframe_gene_overall(gene,tumor):
    df_ensg= pd.read_csv(gene_name_ENSG,sep='\t')
    result_index = df_ensg[(df_ensg['gene_id_version'].str.lower() == gene.lower()) | (df_ensg['gene_id'].str.lower() == gene.lower()) | (df_ensg['gene_symbol'].str.lower() == gene.lower())].index
    #result_index = df_ensg[(df_ensg['gene_id_version'] == gene) | (df_ensg['gene_id'] == gene) | (df_ensg['gene_symbol'] == gene)].index
    if not result_index.empty:
        gene_version=df_ensg.loc[result_index[0],'gene_id_version']
        indice=int(result_index[0])
        gene_dataframe_FPKM_tumor=os.path.join(gene_dataframe_FPKM,"Dataframe_FPKM_"+tumor+".csv")
        df=pd.read_csv(gene_dataframe_FPKM_tumor)
        df=df.set_index("gene_id")
        return(df,gene_version)
    
    
    with open(miRNA_name) as f:
        miRNA_names = f.read().splitlines()
        for mirna in miRNA_names:
            if mirna.lower() == gene.lower():
                df=pd.read_csv(miRNA_dataframe)
                df=df.set_index('miRNA_ID')
                return(df,mirna)
                    
    if gene in open(protein_name).read().split("\n"):
        df=pd.read_csv(protein_dataframe)
        df=df.set_index('peptide_target')
        return (df,gene)
    

    with open(protein_name) as f:
        protein_names = f.read().splitlines()
        for protein in protein_names:
            if protein.lower() == gene.lower():  # Use case-insensitive comparison to find if gene is in the list
                df=pd.read_csv(protein_dataframe)
                df=df.set_index('peptide_target')
                return (df,protein)
        
    print("analysis is not available for the entered name")
    return 0



def dataframe_OStime(tumor,column):
    dfclinic=pd.read_csv(clinical_OS)  
    OS=(dfclinic[['bcr_patient_barcode',column]]) 
    df1_mask=dfclinic['type']==tumor
    OS=dfclinic[df1_mask]
    OS=(OS[['bcr_patient_barcode',column]]) 
    OS=OS.set_index('bcr_patient_barcode')
    #display(OS)
    OS=OS.dropna()
    return(OS)




def overall_survival_analysis(m,tumor,cartella,df1,OS1,gene,gene_input):
    
    
    i1=df1.loc[m,:] > df1.loc[m,:].median()
    i2 = df1.loc[m,:] < df1.loc[m,:].median() 
    
    kmf = KaplanMeierFitter()
    
    #results = logrank_test((OS1[i1]), (OS1[i2]),list(df1.loc[m,i1]),list(df1.loc[m,i2]), alpha=.95)
    results = logrank_test(OS1[i1], OS1[i2], event_observed_A=df1.loc[m, i1], event_observed_B=df1.loc[m, i2])

    #print(results.p_value )
    if results.p_value < 1:       
        print(results.p_value)
        os.makedirs(cartella)
       
        # Creiamo il grafico Kaplan-Meier con miglioramenti
        kmf.fit(OS1[i1], list(df1.loc[m, i1]), label="Higher expression")
        a1 = kmf.plot(color='red', linestyle='-', linewidth=1,ci_show=False)  
        kmf.fit(OS1[i2], list(df1.loc[m, i2]), label="Lower expression")
        kmf.plot(ax=a1, color='blue', linestyle='--', linewidth=1,ci_show=False)  

        # Aggiungiamo una griglia
        plt.grid(True, which='both', linestyle='--', linewidth=0.5)
        # Rimuovere i bordi del grafico
        a1.spines['top'].set_visible(False)   # Rimuove il bordo superiore
        a1.spines['right'].set_visible(False) # Rimuove il bordo destro
        

        # Aggiungiamo etichette e titolo
        plt.xlabel('Time (days)', fontsize=10)
        plt.ylabel('Survival Probability', fontsize=10)
        plt.title('Kaplan-Meier Survival Curves', fontsize=12)

        # Posizioniamo la legenda in alto a destra, senza bordo
        plt.legend(loc='best', frameon=False)

      
     
        #print(cartella+"/overallsurvival_"+gene+"_"+tumor+".jpeg")
        plt.savefig(cartella+"/overallsurvival_"+gene_input+"_"+tumor+".jpeg",dpi=300)
        
    else:
       
        print("Not enough survival data to calculate the analysis")
        return(0)
  



def overall_survival_analysis_pathway(m,tumor,df1,OS1,cartella):
    # if m not in df1.index:
    #     return('search is not available for the name entered')
    # else:
    i1=df1.loc[m,:] > df1.loc[m,:].median()
    i2 = df1.loc[m,:] < df1.loc[m,:].median() 
    
    kmf = KaplanMeierFitter()
    

    
    #results = logrank_test((OS1[i1]), (OS1[i2]),list(df1.loc[m,i1]),list(df1.loc[m,i2]), alpha=.95)
    results = logrank_test(OS1[i1], OS1[i2], event_observed_A=df1.loc[m, i1], event_observed_B=df1.loc[m, i2])

    
    if results.p_value < 1:
        os.mkdir(cartella)
        print(results.p_value)
    

        kmf.fit(OS1[i1], list(df1.loc[m, i1]), label="Higher expression")
        a1 = kmf.plot(color='red', linestyle='-', linewidth=1,ci_show=False)  
        kmf.fit(OS1[i2], list(df1.loc[m, i2]), label="Lower expression")
        kmf.plot(ax=a1, color='blue', linestyle='--', linewidth=1,ci_show=False)  

        # Aggiungiamo una griglia
        plt.grid(True, which='both', linestyle='--', linewidth=0.5)
        # Rimuovere i bordi del grafico
        a1.spines['top'].set_visible(False)   # Rimuove il bordo superiore
        a1.spines['right'].set_visible(False) # Rimuove il bordo destro
        

        # Aggiungiamo etichette e titolo
        plt.xlabel('Time (days)', fontsize=10)
        plt.ylabel('Survival Probability', fontsize=10)
        plt.title('Kaplan-Meier Survival Curves', fontsize=12)

        # Posizioniamo la legenda in alto a destra, senza bordo
        plt.legend(loc='best', frameon=False)

        
        plt.savefig(cartella+"/"+m+"_"+tumor+".png",dpi=300)

    else:
        print("Not enough survival data to calculate the analysis")
        return(0)




def open_gsva_df(tumor):
    file=os.path.join(os_pathway,'gsva_'+tumor+'.csv')
    df=pd.read_csv(file)
    df=df.set_index('Unnamed: 0')
    df.columns= [x.replace(".","-") for x in df.columns]
    return(df)


#######################################################################################

#Deseq2_mutation

def mapping_ensg():
    return(mapped_ENSG)