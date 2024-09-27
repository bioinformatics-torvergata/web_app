import shutil
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import ranksums
import os
import sys
from function_new_prova import copyfile
import numpy as np


#################################    
#             MAIN              #
#################################   

if __name__ == "__main__":
    tumor= sys.argv[1]
    dir= sys.argv[2]
    dir_saveresults=sys.argv[3]

    
    #copia i file nella cartella dei risulati dir_saveresults
    copyfile(tumor,dir,dir_saveresults)