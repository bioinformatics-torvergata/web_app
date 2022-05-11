import pandas as pd
import sys

def converti_name_gene(gene):
    df=pd.read_csv("/mnt/data/notturno/overall_survival_rnainter/gene_symbol_filtrato.csv")
    
    if gene[:4]!='ENSG':
        return(gene) #gene_symbol      
    else:
        df=df.set_index("ENSG")
        name=(df.loc[gene])
        return name['gene_symbol']


def controllo_interazione(gene):
    #apro df interazioni strong mirna-mrna di rnainter
    df=pd.read_csv("/mnt/data/notturno/overall_survival_rnainter/rnainter_homosapiens_strong.csv")
    
    if gene in list(df['Interactor2.Symbol']):
        df= df.set_index('Interactor2.Symbol')
        df1=(df.loc[gene])
        lista=(list(df1['Interactor1.Symbol']))
        
        lista= [x for x in lista if str(x) != 'nan']
        
        return(lista)
    else: 
        print('questo gene non è presente nel file di interazione')

  

gene= sys.argv[1]
gene=converti_name_gene(gene)
lista_interattori=controllo_interazione(gene)
stringa=','.join(lista_interattori)
print(stringa)
