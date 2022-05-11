import pandas as pd
import sys
import os

tumor= sys.argv[1]


path='/mnt/data/notturno/web_app/webserver/rolls/static/media/deseq2/'
dir=os.listdir(path)

listafeature=[]
for feature in dir:
    #print(feature)
    pathfeat=path+(feature)
    listatumor=os.listdir(pathfeat)
    #print(feature, listatumor)
    if tumor in (listatumor):
        listafeature.append(feature)

stringa=','.join(listafeature)
print(stringa)
      
    