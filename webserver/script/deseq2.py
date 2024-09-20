import shutil
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import ranksums
import os
import sys
from function_new_prova import copyfile,plotly_volcano


       
#################################    
#             MAIN              #
#################################   

if __name__ == "__main__":
    tumor= sys.argv[1]
    dir= sys.argv[2]
    dir_saveresults=sys.argv[3]

    
    #copia i file nella cartella dei risulati dir_saveresults
    copyfile(tumor,dir,dir_saveresults)

    #generare file html plotly per volcano plot nella cartella risultati
    path_result=os.path.join(dir_saveresults,'result_'+tumor+'.txt')
    print(path_result)
    df=pd.read_csv(path_result,sep="\t")

    #plt
    plotly_volcano(df,dir_saveresults,tumor)

   


