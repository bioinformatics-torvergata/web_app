import shutil
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import ranksums
import os
import sys
from function_new_prova import copyfile,plotly_volcano
import numpy as np

       
#################################    
#             MAIN              #
#################################   

if __name__ == "__main__":
    tumor= sys.argv[1]
    dir= sys.argv[2] #path file result.txt
    cartella=sys.argv[3] #where save html
    

    df=pd.read_csv(dir,sep="\t")
    df['padj'] = pd.to_numeric(df['padj'], errors='coerce')
    df['padj']=np.log10(df['padj'])*(-1)
    #plt
    plotly_volcano(df,cartella,tumor)

   


