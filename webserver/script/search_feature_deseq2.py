import pandas as pd
import sys
import os
import configparser
from pathlib import Path


tumor= sys.argv[1]

config = configparser.ConfigParser()
script_dir = Path(__file__).parent

# Costruisci il percorso relativo al file di configurazione
config_file = script_dir.parent.parent / 'webserver' / 'webserver' / 'conf.ini'

config.read(config_file) 
#directory base
output_data = config['Paths']['output_data']
base_dir = config.get('Paths', 'base_dir', fallback='')

path_deseq=os.path.join(base_dir, 'deseq2')

dir = [f for f in os.listdir(path_deseq) if not f.startswith('.')]


listafeature=[]
for feature in dir:
  
    pathfeat=path_deseq+'/'+str(feature)
    listatumor=os.listdir(pathfeat)
    if tumor in (listatumor):
        listafeature.append(feature)

stringa=','.join(listafeature)
print(stringa)
      
    