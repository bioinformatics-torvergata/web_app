#!/usr/bin/env python
# coding: utf-8

import pandas as pd
from decimal import Decimal
import sys
import os
import numpy as np
from function_new_prova import plotly_volcano,mapping_ensg

if __name__ == "__main__":
    tumor= sys.argv[1]
    file_result= sys.argv[2] #"/mnt/data/notturno/deseq2/tobacco_smoking_history/BLCA/result_BLCA.txt"
    dir_saveresults=sys.argv[3]
    
    #generare file html plotly per volcano plot nella cartella risultati
    path_result=os.path.join(dir_saveresults,file_result)
    print(path_result)
    df=pd.read_csv(path_result,sep="\t",header=None,skiprows=1)

    df.columns=["ENSG","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"]

    file_mapped=mapping_ensg()
    ens2sym=pd.read_csv(file_mapped,sep="\t",index_col=0)


    GS=[]
    for i in df.ENSG:
        if i in ens2sym.index:
            GS.append(ens2sym.loc[i,"gene_symbol"])
        else:
            GS.append(i)


    df["GeneSymbol"]=GS
    col=["GeneSymbol","ENSG","log2FoldChange","padj"]

    df=df.loc[:,col].set_index("ENSG")

    #df.padj=['%.2E' % Decimal(x) for x in df.padj]
    df.log2FoldChange=[round(x, 3) for x in df.log2FoldChange]

    
   

    df['padj']=np.log10(df['padj'])*(-1)
    #plt
    plotly_volcano(df,dir_saveresults,tumor)

    #sovrascrive il file di risultati
    df.padj=['%.2E' % Decimal(x) for x in df.padj]
    df.to_csv(path_result,sep="\t")



