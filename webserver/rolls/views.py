from django.shortcuts import render
from django.http import HttpResponse
from subprocess import run,PIPE
import sys
from matplotlib import image
from rolls.forms import Gene, deconvolution_form,FormMutationChoice,featuremutationform,formcorrelation,Analisiform_protein,Analisiformcompleto_protein,formSurvival,Analisiform, Deseq2form, Analisiform1, Analisiformcompleto, Analisi_interaction,Analisipath,tumorGeneform,FormTumorMutation
import os
from django.http import StreamingHttpResponse
from wsgiref.util import FileWrapper
import mimetypes
from shutil import make_archive
import time
import os.path
from django.conf import settings
import shutil
import csv
import pandas as pd
from decimal import Decimal
from django.http import JsonResponse
from rolls.models import Gene,Pathway,Protein, Gene_symbol
from django.http import FileResponse, HttpResponseNotFound
import subprocess


import json
import configparser
from pathlib import Path

#richiamo conf.ini
config = configparser.ConfigParser()
script_dir = Path(__file__).parent

# Costruisci il percorso relativo al file di configurazione
config_file = script_dir.parent.parent / 'webserver' / 'webserver' / 'conf.ini'

config.read(config_file) 
#directory base
output_data = config['Paths']['output_data']
base_dir = config.get('Paths', 'base_dir', fallback='')
output_data_Table=config['Paths']['output_data_Table']

def get_full_path(relative_path):
    return os.path.join(base_dir, relative_path)

parametri={
    'patient_status':'Tumor vs Ctrl',
    'gender':'Female vs Male',
    'pathologic_stage': 'Stage III-IV vs Stage I-II',
    'radiation_therapy': 'YES vs NO',
    'diabetes': 'YES vs NO',
    'tobacco_smoking_history': 'Smoker vs Non_Smoker',
    'menopause_status':'Post-menopause vs Pre-menopause',
    'alcohol_history_documented': 'YES vs NO',
    'age_at_initial_pathologic_diagnosis': 'Above the median vs Below the median' ,
}




def rolls(request):
    count = get_counter()
    return render(request, 'rolls/home.html',{'count': count,})

def documentation(request):
    count = get_counter()
    return render(request, 'rolls/documentation.html',{'count': count,})

def read_table(file_path):
    txt_data = []

    df = pd.read_csv(file_path, sep='\t',dtype=str)
    #print(df)

    txt_data.append(df.columns.tolist())
    
    # Itera sulle righe del DataFrame e aggiungi ogni riga come lista
    for index, row in df.iterrows():
        txt_data.append(row.tolist())
    return(txt_data)


def dataset(request):
    count = get_counter()
    file_path = os.path.join(output_data_Table,'table','campioni_TCGA.txt')
    file_path_genetype=os.path.join(output_data_Table,'table','gene_type.txt')
    file_feature_tumor=os.path.join(output_data_Table,'table','Features_tumor.txt')
    file_tumor=os.path.join(output_data_Table,'table','tumor_abbreviations.txt')
    file_parameters=os.path.join(output_data_Table,'table','feature_parameters.txt')

    txt_data_clinical=read_table(file_path)

    txt_data_typegene=read_table(file_path_genetype)   
    txt_feature_tumor= read_table(file_feature_tumor)
    txt_tumor= read_table(file_tumor) 
    txt_parameters=read_table(file_parameters)
    return render(request, 'rolls/dataset.html', {
        'dati': txt_data_clinical, 
        'dati_genetype':txt_data_typegene,
        'dati_feature':txt_feature_tumor,
        'dati_tumor':txt_tumor,
        'dati_parametri':txt_parameters,
        'count': count,
 })


def contact(request):
    count = get_counter()
    return render(request, 'rolls/contact.html',{'count': count,})

def yourdataset(request):
    return render(request,'rolls/yourdataset.html')


##########################################
# Function for autocomplete
def gene_suggestions(request):
    if 'term' in request.GET:
        qs = Gene.objects.filter(gene__istartswith=request.GET.get('term'))
        genes = sorted(qs.values_list('gene', flat=True), key=len)[:10]  # Ordina per lunghezza e limita a 10
        return JsonResponse(genes, safe=False)
    
    
def pathway_suggestions(request):
    if 'term' in request.GET:
        qs = Pathway.objects.filter(pathway__icontains=request.GET.get('term'))[:20]  # Limita i risultati a 10
        pathways = list(qs.values_list('pathway', flat=True))
        return JsonResponse(pathways, safe=False)
    
def protein_suggestions(request):
    if 'term' in request.GET:
        qs = Protein.objects.filter(protein__istartswith=request.GET.get('term'))
        protein = sorted(qs.values_list('protein', flat=True), key=len)[:10]  # Ordina per lunghezza e limita a 20
        return JsonResponse(protein, safe=False)
    
def gene_symbol_suggestions(request):
    if 'term' in request.GET:
        qs = Gene_symbol.objects.filter(gene_symbol__istartswith=request.GET.get('term'))
        gene_symbol = sorted(qs.values_list('gene_symbol', flat=True), key=len)[:10]  # Ordina per lunghezza e limita a 20
        return JsonResponse(gene_symbol, safe=False)
##########################################



#contatore per visualizzare quante analisi sono state effettuate
counter_file = 'analysis_counter.txt'
def get_counter():
    # Se il file non esiste, ritorna 0 come valore di default
    if not os.path.exists(counter_file):
        return 0
    
    # Leggi il valore del contatore dal file
    with open(counter_file, 'r') as f:
        count = int(f.read().strip())
    
    return count

# Funzione per aggiornare il contatore
def increment_counter():
        # Se il file non esiste, crealo con valore iniziale 0
        if not os.path.exists(counter_file):
            with open(counter_file, 'w') as f:
                f.write('0')
        
        # Leggi il valore corrente del contatore, incrementa e riscrivilo
        with open(counter_file, 'r+') as f:
            count = int(f.read().strip())
            count += 1
            f.seek(0)
            f.write(str(count))
            f.truncate()
        return count
####################################################################################



TUMOR_FEATURE_MAPPING = {
    "ACC":['age_at_initial_pathologic_diagnosis','gender','pathologic_stage','radiation_therapy'],
    "BLCA":['age_at_initial_pathologic_diagnosis','gender','pathologic_stage','patient_status','radiation_therapy','tobacco_smoking_history'],
    "BRCA":['age_at_initial_pathologic_diagnosis','menopause_status','pathologic_stage','patient_status','radiation_therapy'],
    "CESC":['age_at_initial_pathologic_diagnosis','menopause_status','radiation_therapy','tobacco_smoking_history'],
    "CHOL":['age_at_initial_pathologic_diagnosis','gender','pathologic_stage'],
    "COAD":['age_at_initial_pathologic_diagnosis','gender','pathologic_stage','patient_status','radiation_therapy'],
    "DLBC":['age_at_initial_pathologic_diagnosis','gender','radiation_therapy'],
    "ESCA":['age_at_initial_pathologic_diagnosis','alcohol_history_documented','gender','pathologic_stage','radiation_therapy','tobacco_smoking_history'],
    "GBM":['age_at_initial_pathologic_diagnosis','gender','radiation_therapy'],
    "HNSC":['age_at_initial_pathologic_diagnosis','alcohol_history_documented','gender','pathologic_stage','patient_status','radiation_therapy','tobacco_smoking_history'],
    "KICH":['age_at_initial_pathologic_diagnosis','gender','pathologic_stage','patient_status','tobacco_smoking_history'],
    "KIRC":['age_at_initial_pathologic_diagnosis','gender','pathologic_stage','patient_status','radiation_therapy','tobacco_smoking_history'],
    "KIRP":['age_at_initial_pathologic_diagnosis','gender','pathologic_stage','patient_status','radiation_therapy','tobacco_smoking_history'],
    "LGG":['age_at_initial_pathologic_diagnosis','gender','radiation_therapy'],
    "LIHC":['age_at_initial_pathologic_diagnosis','gender','pathologic_stage','patient_status','radiation_therapy'],
    "LUAD":['age_at_initial_pathologic_diagnosis','gender','pathologic_stage','patient_status','radiation_therapy','tobacco_smoking_history'],
    "LUSC":['age_at_initial_pathologic_diagnosis','gender','pathologic_stage','patient_status','radiation_therapy','tobacco_smoking_history'],
    "MESO":['age_at_initial_pathologic_diagnosis','gender','pathologic_stage','radiation_therapy'],
    "OV":['age_at_initial_pathologic_diagnosis','radiation_therapy'],
    "PAAD":['age_at_initial_pathologic_diagnosis','alcohol_history_documented','gender','pathologic_stage','radiation_therapy','tobacco_smoking_history'],
    "PCPG":['age_at_initial_pathologic_diagnosis','gender','radiation_therapy'],
    "PRAD":['age_at_initial_pathologic_diagnosis','patient_status','radiation_therapy'],
    "READ":['age_at_initial_pathologic_diagnosis','gender','pathologic_stage','radiation_therapy'],
    "SARC":['age_at_initial_pathologic_diagnosis','gender','radiation_therapy'],
    "SKCM":['age_at_initial_pathologic_diagnosis','gender','pathologic_stage','radiation_therapy'],
    "STAD":['age_at_initial_pathologic_diagnosis','gender','pathologic_stage','patient_status','radiation_therapy'],
    "TGCT":['age_at_initial_pathologic_diagnosis','pathologic_stage','radiation_therapy'],
    "THCA":['age_at_initial_pathologic_diagnosis','gender','pathologic_stage','patient_status','radiation_therapy'],
    "THYM":['age_at_initial_pathologic_diagnosis','gender','radiation_therapy'],
    "UCEC":['age_at_initial_pathologic_diagnosis','diabetes','menopause_status','patient_status','radiation_therapy'],
    "UCS":['age_at_initial_pathologic_diagnosis','diabetes','radiation_therapy'],
    "UVM":['age_at_initial_pathologic_diagnosis','gender','pathologic_stage','radiation_therapy'],
}


def get_features(request):
    # Ottieni il valore del tumore selezionato
    tumor = request.GET.get('tumor')
    
    # Ottieni le feature corrispondenti dal dizionario
    features = TUMOR_FEATURE_MAPPING.get(tumor, [])
    
    # Restituisci le feature come tuple (valore, etichetta)
    feature_list = [(feature, feature) for feature in features]
    
    return JsonResponse(feature_list, safe=False)



def download_zip(request, subdirectory, zip_name):
    # Componi il percorso completo del file ZIP usando il terzo parametro
    zip_file_path = os.path.join(output_data, subdirectory, zip_name)
    print(zip_file_path)

    # Verifica se il file esiste e se è all'interno di una directory sicura
    # Controlla se il file è all'interno della directory 'output_data' per prevenire attacchi di directory traversal
    if os.path.commonprefix([os.path.realpath(zip_file_path), os.path.realpath(output_data)]) == os.path.realpath(output_data):
        if os.path.exists(zip_file_path):
            return FileResponse(open(zip_file_path, 'rb'), as_attachment=True, filename=zip_name)
        else:
            return HttpResponseNotFound("File non trovato.")
    else:
        return HttpResponseNotFound("Accesso non autorizzato al file.")


############################################
#                                          #
#             TRANSCRIPTOMIC               #
#                                          #
############################################

# analisi: 
#1. Deseq2 (OK)
#2. Differential expression single tumor (OK)
#3. Differential expression all tumor (OK)


#############      DESEQ2 analisi     ###############  

def choseimage(tumor,pathfiles,dir_saveresults):
    print(' i file li pesca da qui ', dir_saveresults)
    files=os.listdir(dir_saveresults)
    print(files)
    filelist=[]
    for file in files:
        print(file)
        if tumor in file:
            if 'jpg' in file:
                if 'EnhancedVolcano' in file:
                    enhancedimage=file
                if 'heatmap' in file:
                    heatmap=file
                if 'PCA' in file: 
                    pca=file
                if 'TopGeni' in file:
                    topgeni=file
                else:
                    filelist.append(file)
            
            if 'html' in file:
                image_plotly=file


    if len(filelist)>0:
        print(filelist)
        return(enhancedimage,heatmap, pca,topgeni,image_plotly)
    else:
        return()


def read_table_deseq(file_path):
    txt_data = []
    
    df = pd.read_csv(file_path, sep='\t') 
    df['padj'] = pd.to_numeric(df['padj'], errors='coerce')
    #df = df.sort_values(by='padj', ascending=True).head(500)
    df=df[df['padj'] < 0.05].sort_values(by='padj', ascending=True) #.head(500)

    df.padj=['%.2E' % Decimal(x) for x in df.padj]
    
    # Aggiungi l'intestazione (nomi delle colonne, incluso l'indice) alla lista
    txt_data.append(df.columns.tolist())
    
    # Itera sulle righe del DataFrame e aggiungi ogni riga come lista
    for index, row in df.iterrows():
        txt_data.append(row.tolist())
    return(txt_data)


def deseq2(request):
    count = get_counter()
    if request.method == 'POST':

            form = Deseq2form(request.POST)
            tumor=request.POST['tumor'] 
            feature=request.POST.get('feature', False)
            dir= os.path.join(base_dir,'deseq2', feature,tumor)
            
            
            if os.path.isdir(dir): 
                inp3=(time.strftime("%Y-%m-%d-%H-%M-%S"))
                dir_saveresults= os.path.join(output_data, inp3)
                print(dir_saveresults)
                os.makedirs(dir_saveresults)

                result_file='result_' + tumor + '_2.txt'
                out=run([sys.executable,'script/deseq2.py',tumor,dir,dir_saveresults,result_file],shell=False, stdout=PIPE)
                count = increment_counter()
                images=choseimage(tumor,dir,dir_saveresults)
                
                file_txt=os.path.join(output_data,inp3,result_file)
                print(file_txt)
                result=read_table_deseq(file_txt)
        


                form=Deseq2form()
                return render(request, 'rolls/deseq2.html', {'form':form, 
                    'feature': feature,
                    'tumor':tumor,
                    'enhancedimage': os.path.join('media/saveanalisi',inp3,images[0]),
                    'images1': os.path.join('media/saveanalisi',inp3,images[1]),
                    'images2': os.path.join('media/saveanalisi',inp3,images[2]),
                    'images3': os.path.join('media/saveanalisi',inp3,images[3]),
                    'image_plotly':os.path.join('media/saveanalisi',inp3,images[4]),
                    'go':'Valid',
                    'parametri': parametri[feature],
                    'table':'media/saveanalisi/'+inp3+'/'+result_file,
                    'dir':inp3,
                    'dati':result,
                    'count': count,
                    })

            else:
                form=Deseq2form()
                return render(request, 'rolls/deseq2.html', {'form':form,
                'feature': feature,
                'tumor':tumor, 
                'count': count,
                'go':'error'})



    form = Deseq2form()       
    return render(request, 'rolls/deseq2.html', {'form':form,'count': count,})


#############      DIFFERENTIAL EXPRESSION SINGLE TUMOR  TRASCRITTOMIC   ############
def diff_exp_single_tumor(request):
    count = get_counter()
    if request.method == 'POST':
            form = Analisiformcompleto(request.POST)
            if form.is_valid():
                gene=form.cleaned_data['gene']
                tumor=form.cleaned_data['tumor']
                feature=form.cleaned_data['feature']
                
                time_dir=(time.strftime("%Y-%m-%d-%H-%M-%S"))
                dir= os.path.join(output_data, time_dir)
                control=''
                os.makedirs(dir)
                out=run([sys.executable,'script/Differential_expression_boxplot_plotly_new.py',gene,tumor,feature,dir,control],shell=False, stdout=PIPE)
                debug_error=out.stdout.decode().strip()
                print((debug_error))

                count = increment_counter()
                
                if debug_error=='0':
                    form=Analisiformcompleto()
                    return render(request, 'rolls/diff_exp_single_tumor.html', {'form':form,
                    'formresult': out.stdout.decode('ascii'),
                    'gene':gene,
                    'count': count,
                    'go':'error_name'})
                
                if debug_error=='2':
                    form=Analisiformcompleto()
                    return render(request, 'rolls/diff_exp_single_tumor.html', {'form':form,
                    'formresult': out.stdout.decode('ascii'),
                    'gene':gene,
                    'tumor':tumor,
                    'feature':feature,
                    'count': count,
                    'go':'error'})
                else:
                    
                    if os.path.isdir(dir): 
                        files=os.listdir(dir)
                        n=0
                        for file in files:
                            
                            if 'html' in file :
                                image=os.path.join('media/saveanalisi',time_dir,file)
                                n+=1

                            if 'Result' in file:
                               
                                result_data=os.path.join(output_data,dir,file)
                                result=read_table_comma(result_data)
                                n+=1
                            # image='/media/saveanalisi/'+time_dir+'/'+file
                        
                        if n>0: 
                            #zip folder analisi -> results.zip
                            folder_to_zip=os.path.join(dir,"results.zip")
                            subprocess.run(["zip", "-r", folder_to_zip, "."],cwd=dir) 

                            form=Analisiformcompleto()
                            return render(request, 'rolls/diff_exp_single_tumor.html', {
                                        'form':form, 
                                        'formresult': out.stdout.decode('ascii'),
                                        'image': image,
                                        'go':'Valid',
                                        'gene':gene,
                                        'tumor':tumor,
                                        'feature':feature,
                                        'parametri': parametri[feature],
                                        'dati':result,
                                        'count': count,
                                        'dir':time_dir,
                                        'parametri': parametri[feature],
                                        })
                        if n==0:
                        
                            form=Analisiformcompleto()
                            return render(request, 'rolls/diff_exp_single_tumor.html', {'form':form,
                            'formresult': out,
                            'feature': feature,
                            'tumor':tumor, 
                            'gene':gene,
                            'count': count,
                            'go':'error'})
                


    form=Analisiformcompleto()
    return render(request, 'rolls/diff_exp_single_tumor.html', {'form':form, 'count': count,})


#############      DIFFERENTIAL EXPRESSION ANALYSIS ALL TUMOR FOR FEATURE TRASCRITTOMIC   #############
def read_table_all_tumor(file_path):
    txt_data = []

    df = pd.read_csv(file_path, sep='\t',dtype=str) 

    df['p-value'] = pd.to_numeric(df['p-value'], errors='coerce')
    df = df.dropna(subset=['p-value'])
    df = df.sort_values(by='p-value', ascending=True)
    
    df['p-value']=['%.2E' % Decimal(x) for x in df['p-value']]
     # Aggiungi l'intestazione (nomi delle colonne, incluso l'indice) alla lista
    txt_data.append(df.columns.tolist())
    
    # Itera sulle righe del DataFrame e aggiungi ogni riga come lista
    for index, row in df.iterrows():
        txt_data.append(row.tolist())
    return(txt_data)
    


def differential_expression(request):
    count = get_counter()
    if request.method == 'POST':
        form = Analisiform(request.POST)
        if form.is_valid():
            gene=form.cleaned_data['gene']
            feature=form.cleaned_data['feature']
            
            inp3=(time.strftime("%Y-%m-%d-%H-%M-%S"))
            control=''
            dir= os.path.join(output_data, inp3)
            out=run([sys.executable,'script/boxplot_all_tumor_giusto_new.py',gene,feature,dir,control],shell=False, stdout=PIPE)
            print(out)
            count = increment_counter()
            debug_error=out.stdout.decode().strip()
            print((debug_error))
            
            if debug_error=='0':
                form=Analisiform()
                return render(request, 'rolls/differential_expression.html', {'form':form,
                'formresult': out.stdout.decode('ascii'),
                'gene':gene,
                'count': count,
                'go':'error_name'})
            
            if os.path.isdir(dir): 
                files=os.listdir(dir)
                for file in files:
                    if 'html' in file:
                        image=os.path.join('media/saveanalisi',inp3,file)
                        #image='media/saveanalisi/'+inp3+'/'+file

                    if 'result' in file:
                        result_data=os.path.join(output_data,inp3,file)
                        result=read_table_all_tumor(result_data)
                
                #zip folder analisi -> results.zip
                folder_to_zip=os.path.join(dir,"results.zip")
                subprocess.run(["zip", "-r", folder_to_zip, "."],cwd=dir)   

                form=Analisiform()
                return render(request, 'rolls/differential_expression.html', {
                    'form':form, 
                    'formresult': out.stdout.decode('ascii'),
                    'image': image,
                    'go':'Valid',
                    'gene':gene,
                    'feature':feature,
                    'parametri': parametri[feature],
                    #'dir':'media/saveanalisi/'+inp3+'/result.txt', #download singola tabella
                    'dir':inp3,
                    'dati':result,
                    'count': count,
                    })
            else:
                form=Analisiform()
                return render(request, 'rolls/differential_expression.html', {
                    'form':form,
                    'formresult': out.stdout.decode('ascii'),
                    'feature': feature,
                    'gene':gene,
                    'count': count, 
                    'go':'error'})


    form=Analisiform()
    return render(request, 'rolls/differential_expression.html', {'form':form,'count': count,})





############################################
#                                          #
#                PROTEOMIC                 #
#                                          #
############################################

# analisi: 
#1. Differential expression single tumor (OK)
#2. Differential expression all tumor (OK -controllare bug NAN)

############     DIFFERENTIAL EXPRESSION SINGLE TUMOR PROTEOMIC   ############
def diff_exp_single_tumor_protein(request):
    count = get_counter()
    if request.method == 'POST':
            form = Analisiformcompleto_protein(request.POST)
            if form.is_valid():
                gene=form.cleaned_data['protein']
                tumor=form.cleaned_data['tumor']
                feature=form.cleaned_data['feature']
                control='protein'
                time_dir=(time.strftime("%Y-%m-%d-%H-%M-%S"))
                dir= os.path.join(output_data, time_dir)
               # dir=os.path.join(settings.BASE_DIR, 'rolls', 'static', 'media', 'saveanalisi', time_dir)
                os.makedirs(dir)
                out=run([sys.executable,'script/Differential_expression_boxplot_plotly_new.py',gene,tumor,feature,dir,control],shell=False, stdout=PIPE)
                print(out)
                count = increment_counter()
                debug_error=out.stdout.decode().strip()
                print((debug_error))
                
                if debug_error=='0':
                    form=Analisiformcompleto()
                    return render(request, 'rolls/diff_exp_single_tumor_protein.html', {'form':form,
                    'formresult': out.stdout.decode('ascii'),
                    'gene':gene,
                    'count': count,
                    'go':'error_name'})
                
                if debug_error=='2':
                    form=Analisiformcompleto()
                    return render(request, 'rolls/diff_exp_single_tumor_protein.html', {'form':form,
                    'formresult': out.stdout.decode('ascii'),
                    'gene':gene,
                    'tumor':tumor,
                    'feature':feature,
                    'count': count,
                    'go':'error'})
                else:
                
                    if os.path.isdir(dir): 
                        files=os.listdir(dir)
                        n=0
                        for file in files:
                            
                            if 'html' in file :
                                image=os.path.join('media/saveanalisi',time_dir,file)
                                n+=1

                            if 'Result' in file:
                               
                                result_data=os.path.join(output_data,dir,file)
                                result=read_table_comma(result_data)
                                n+=1
                           
                        if n>0: 
                            #zip folder analisi -> results.zip
                            folder_to_zip=os.path.join(dir,"results.zip")
                            subprocess.run(["zip", "-r", folder_to_zip, "."],cwd=dir)

                            form=Analisiformcompleto_protein()
                            return render(request, 'rolls/diff_exp_single_tumor_protein.html', {
                                'form':form, 
                                'formresult': out.stdout.decode('ascii'),
                                'image': image,
                                'go':'Valid',
                                'gene':gene,
                                'tumor':tumor,
                                'feature':feature,
                                'parametri': parametri[feature],
                                'count': count,
                                'dir': time_dir,
                                'dati':result,})
                    else:
                        form=Analisiformcompleto_protein()
                        return render(request, 'rolls/diff_exp_single_tumor_protein.html', {'form':form,
                        'feature': feature,
                        'tumor':tumor, 
                        'gene':gene,
                        'count': count,
                        'go':'error'})


    form=Analisiformcompleto_protein()
    return render(request, 'rolls/diff_exp_single_tumor_protein.html', {'form':form,'count': count,})

#############    DIFFERENTIAL EXPRESSION ANALYSIS ALL TUMOR FOR FEATURE   PROTEOMIC   ############# 
def differential_expression_protein(request):
    count = get_counter()
    if request.method == 'POST':
        form = Analisiform_protein(request.POST)
        if form.is_valid():
            gene=form.cleaned_data['protein']
            feature=form.cleaned_data['feature']
            
            inp3=(time.strftime("%Y-%m-%d-%H-%M-%S"))
            control='protein'
            dir= os.path.join(output_data, inp3)
            out=run([sys.executable,'script/boxplot_all_tumor_giusto_new.py',gene,feature,dir,control],shell=False, stdout=PIPE)
            print(out)
            count = increment_counter()
            debug_error=out.stdout.decode().strip()
            print(debug_error)
            
            if debug_error=='0':
                form=Analisiform_protein()
                return render(request, 'rolls/differential_expression_protein.html', {'form':form,
                'formresult': out.stdout.decode('ascii'),
                'gene':gene,
                'count': count,
                'go':'error_name'})
            
            if os.path.isdir(dir): 
                files=os.listdir(dir)
                for file in files:
                    if 'html' in file:
                        image=os.path.join('media/saveanalisi',inp3,file)
                        #image='media/saveanalisi/'+inp3+'/'+file
                    if 'txt' in file:
                        result_data=os.path.join(output_data,inp3,file)
                        result=read_table_all_tumor(result_data)
                form=Analisiform_protein()

                #zip folder analisi -> results.zip
                folder_to_zip=os.path.join(dir,"results.zip")
                subprocess.run(["zip", "-r", folder_to_zip, "."],cwd=dir) 

                return render(request, 'rolls/differential_expression_protein.html', {
                    'form':form, 
                    'formresult': out.stdout.decode('ascii'),
                    'image': image,
                    'go':'Valid',
                    'gene':gene,
                    'feature':feature,
                    'parametri': parametri[feature],
                    #'dir':'media/saveanalisi/'+inp3+'/result.txt', #download singola tabella
                    'dir':inp3,
                    'dati':result,
                    'count': count,
                    'parametri': parametri[feature],
                    })
            else:
                form=Analisiform_protein()
                return render(request, 'rolls/differential_expression_protein.html', {
                    'form':form,
                    'formresult': out.stdout.decode('ascii'),
                    'feature': feature,
                    'gene':gene, 
                    'count': count,
                    'go':'error'})


    form=Analisiform_protein()
    return render(request, 'rolls/differential_expression_protein.html', {'form':form,'count': count,})



############################################
#                                          #
#                SURVIVAL                  #
#                                          #
############################################

# analisi:
#1. Overall survival (OK - inserire pvalue nel grafico ? numerosità campionaria?)
#2. Overall survival con pathway activity score (OK - inserire pvalue nel grafico? numerosità campionaria?)
#3. Survival with gene mutation status (OK)


########### OVERALL SURVIVAL ###########
def overall_survival(request):
    count = get_counter()
    if request.method == 'POST':
        form = formSurvival(request.POST)
        if form.is_valid():
            gene=form.cleaned_data['gene']
            tumor=form.cleaned_data['tumor']
            
            methods=form.cleaned_data['Methods']

            print(gene, tumor,methods)
            inp3=(time.strftime("%Y-%m-%d-%H-%M-%S"))
            dir= os.path.join(output_data, inp3)
            
            out=run([sys.executable,'script/overall_survival.py',gene,tumor,dir,methods],shell=False, stdout=PIPE)
            print(out)
            count = increment_counter()

            if os.path.isdir(dir): 
                files=os.listdir(dir)
                for file in files:
                    if 'jpeg' in file:
                        image=os.path.join('media/saveanalisi',inp3,file)
                        
                        #zip folder analisi -> results.zip
                        folder_to_zip=os.path.join(dir,"results.zip")
                        subprocess.run(["zip", "-r", folder_to_zip, "."],cwd=dir)

                        form=formSurvival()
                        return render(request, 'rolls/overall_survival.html', {
                            'form':form, 
                            'formresult': out.stdout.decode('ascii'),
                            'image': image,
                            'go':'Valid',
                            'gene':gene ,
                            'tumor':tumor,
                            'method':methods,
                            'count': count,
                            'dir':inp3, 
                            })

            else:
                form=formSurvival()
                return render(request, 'rolls/overall_survival.html', {'form':form,
                'formresult': out.stdout.decode('ascii'),
                'gene':gene,
                'tumor':tumor, 
                'count': count,
                'go':'error'})


    form=formSurvival()
    return render(request, 'rolls/overall_survival.html', {'form':form,'count': count,})


########### Overall survival con pathway activity score ###########
def os_pathway(request):
    count = get_counter()
    if request.method == 'POST':
        form = Analisipath(request.POST)
        if form.is_valid():
            pathway=form.cleaned_data['pathway']
            tumor=form.cleaned_data['tumor']
            method=form.cleaned_data['Methods']
            inp3=(time.strftime("%Y-%m-%d-%H-%M-%S"))
            dir= os.path.join(output_data, inp3)
            out=run([sys.executable,'script/OS_pathway.py',tumor,pathway,dir,method],shell=False, stdout=PIPE)
            print(out)
            count = increment_counter()
            
            if os.path.isdir(dir): 
                files=os.listdir(dir)
                for file in files:
                    if file[-3:]=='png':
                        image='/media/saveanalisi/'+inp3+'/'+file
                
                #zip folder analisi -> results.zip
                folder_to_zip=os.path.join(dir,"results.zip")
                subprocess.run(["zip", "-r", folder_to_zip, "."],cwd=dir) 

                form=Analisipath()
                return render(request, 'rolls/OS_pathway.html', {
                    'form':form, 
                    'formresult': out.stdout.decode('ascii'),
                    'image': image,
                    'go':'Valid',
                    'pathway':pathway ,
                    'tumor':tumor,
                    'count': count,
                    'method':method,
                    'dir':inp3,
                    })
            else:
                form=Analisipath()
                return render(request, 'rolls/OS_pathway.html', {'form':form,
                'pathway':pathway ,
                'tumor':tumor, 
                'count': count,
                'go':'error'})


    form=Analisipath()
    return render(request, 'rolls/OS_pathway.html', {'form':form,'count': count,})


########### Survival with gene mutation status ########### 

def survival_with_gene_mutation_status(request):
    count = get_counter()
    if request.method == 'POST':

        form = tumorGeneform(request.POST)
        if form.is_valid(): 
            tumor=request.POST['tumor'] 
            gene=request.POST['gene'] 
            inp3=(time.strftime("%Y-%m-%d-%H-%M-%S"))
            dir=os.path.join(output_data, inp3)
            # print
            # os.makedirs(dir)
            
            out = subprocess.run(['Rscript', 'script/survival_with_gene_mutation_status.R',tumor,gene,dir], capture_output=True, text=True)
            print(out)
            n=0
            count = increment_counter()
            if os.path.isdir(dir): 
                files=os.listdir(dir)
                for file in files:
                    if 'png' in file:
                        n+=1
                        image=os.path.join('media/saveanalisi',inp3,file)    
                
                        #zip folder analisi -> results.zip
                        folder_to_zip=os.path.join(dir,"results.zip")
                        subprocess.run(["zip", "-r", folder_to_zip, "."],cwd=dir)

                        form=tumorGeneform()
                        return render(request, 'rolls/survival_with_gene_mutation_status.html', {'form':form, 
                            'gene':gene,
                            'tumor':tumor,
                            'image':image,
                            'go':'Valid',
                            'dir':inp3,
                            'count': count,
                            })
                if n==0:
                    form=tumorGeneform()
                    return render(request, 'rolls/survival_with_gene_mutation_status.html', {
                    'form':form,
                    'formresult':'analysis is not available for the entered name',                                                                      
                    'tumor':tumor,
                    'count': count, 
                    'go':'error'})
            else:
                
                form=tumorGeneform()
                return render(request, 'rolls/survival_with_gene_mutation_status.html', {
                'form':form,
                'formresult':'analysis is not available for the entered name',                                                                      
                'tumor':tumor,
                'count': count, 
                'go':'error'})



    form = tumorGeneform()       
    return render(request, 'rolls/survival_with_gene_mutation_status.html', {'form':form,'count': count,})




############################################
#                                          #
#                MUTATION                  #
#                                          #
############################################

#analisi:
#1. tumor_mutation_analysis (OK)
#2. tumor_oncoplot (OK)
#3. somatic_interaction_analysis (OK)
#4. gene mutation analysis (OK)
#5. Differential expression for mutated status - Deseq2 (OK)
#6. Differential mutated gene by clinical feature (OK)


########### TUMOR MUTATION ANALYSES ########### 
def tumor_mutation_analysis(request):
    count = get_counter()
    if request.method == 'POST':

        form = FormTumorMutation(request.POST)
        if form.is_valid(): 
            tumor=request.POST['tumor'] 
            
            inp3=(time.strftime("%Y-%m-%d-%H-%M-%S"))
            dir=os.path.join(output_data, inp3)
            
            os.makedirs(dir)
            
            out = subprocess.run(['Rscript', 'script/tumor_mutation_analysis.R',tumor,dir], capture_output=True, text=True)
            print(out)
            count = increment_counter()
            
            if os.path.isdir(dir): 
                files=os.listdir(dir)
                print(files)
                for file in files:
                    if file[-3:]=='png':
                        if 'summary' in file:
                            image_summary=os.path.join('media/saveanalisi',inp3,file)  
                            
                       
                        if 'Titv' in file:
                            image_titv=os.path.join('media/saveanalisi',inp3,file)
                        

                 #zip folder analisi -> results.zip
                folder_to_zip=os.path.join(dir,"results.zip")
                subprocess.run(["zip", "-r", folder_to_zip, "."],cwd=dir)


                form=FormTumorMutation()
                return render(request, 'rolls/tumor_mutation_analysis.html', {'form':form, 
                    'tumor':tumor,
                    'image_summary':image_summary,
                    'image_titv':image_titv,
                    'go':'Valid',
                    'dir':inp3,
                    'count': count,
                    })

            else:
                form=FormTumorMutation()
                return render(request, 'rolls/tumor_mutation_analysis.html', {'form':form,
                'tumor':tumor, 
                'count': count,
                'go':'error'})



    form = FormTumorMutation()       
    return render(request, 'rolls/tumor_mutation_analysis.html', {'form':form,'count': count,})


########### TUMOR ONCOPLOT  ########### 
def tumor_oncoplot(request):
    count = get_counter()
    if request.method == 'POST':

        form = FormMutationChoice(request.POST)
        if form.is_valid(): 
            tumor=request.POST['tumor'] 
            number=request.POST['number'] 
            inp3=(time.strftime("%Y-%m-%d-%H-%M-%S"))
            dir=os.path.join(output_data, inp3)
            
            os.makedirs(dir)
            
            out = subprocess.run(['Rscript', 'script/tumor_oncoplot.R',tumor,dir,number], capture_output=True, text=True)
            print(out)
            count = increment_counter()
            if os.path.isdir(dir): 
                files=os.listdir(dir)
                print(files)
                for file in files:
                    if file[-3:]=='png':                            
                        if 'oncoplot' in file:
                            image_oncoplot=os.path.join('media/saveanalisi',inp3,file)

                #zip folder analisi -> results.zip
                folder_to_zip=os.path.join(dir,"results.zip")
                subprocess.run(["zip", "-r", folder_to_zip, "."],cwd=dir)

                        
                form=FormMutationChoice()
                return render(request, 'rolls/tumor_oncoplot.html', {'form':form, 
                    'tumor':tumor,
                    'image_oncoplot':image_oncoplot,
                    'go':'Valid',
                    'dir':inp3,
                    'count': count,
                    })

            else:
                form=FormMutationChoice()
                return render(request, 'rolls/tumor_oncoplot.html', {'form':form,
                'tumor':tumor, 
                'count': count,
                'go':'error'})



    form = FormMutationChoice()       
    return render(request, 'rolls/tumor_oncoplot.html', {'form':form,'count': count,})


########### SOMATIC INTERACTION ANALYSIS ########

def read_table_comma(file_path):
    txt_data = []

    df = pd.read_csv(file_path, sep=',',dtype=str) 
    print(df)
    txt_data.append(df.columns.tolist())
    
    # Itera sulle righe del DataFrame e aggiungi ogni riga come lista
    for index, row in df.iterrows():
        txt_data.append(row.tolist())
    return(txt_data)

def somatic_interaction_analysis(request):
    count = get_counter()
    if request.method == 'POST':

        form = FormMutationChoice(request.POST)
        if form.is_valid(): 
            tumor=request.POST['tumor'] 
            number=request.POST['number'] 
            inp3=(time.strftime("%Y-%m-%d-%H-%M-%S"))
            dir=os.path.join(output_data, inp3)
            
            os.makedirs(dir)
            
            out = subprocess.run(['Rscript', 'script/somatic_interaction_analysis.R',tumor,dir,number], capture_output=True, text=True)
            print(out)
            count = increment_counter()
            if os.path.isdir(dir): 
                files=os.listdir(dir)
                print(files)
                for file in files:
                    if file[-3:]=='png':
                        if 'somaticInteractions' in file:
                            image_interact=os.path.join('media/saveanalisi',inp3,file)  
                    if 'results' in file:
                        name=file
                        result_data=os.path.join(output_data,inp3,file)
                        result=read_table_comma(result_data)

                #zip folder analisi -> results.zip
                folder_to_zip=os.path.join(dir,"results.zip")
                subprocess.run(["zip", "-r", folder_to_zip, "."],cwd=dir)

                form=FormMutationChoice()
                return render(request, 'rolls/somatic_interaction_analysis.html', {'form':form, 
                    'tumor':tumor,
                    'number':number,
                    'image_interact':image_interact,
                    'go':'Valid',
                    'dir':inp3,
                    'dati':result,
                    #'dir':'media/saveanalisi/'+inp3+'/'+name, #dowload single table
                    'count': count,
                    })

            else:
                form=FormMutationChoice()
                return render(request, 'rolls/somatic_interaction_analysis.html', {'form':form,
                'tumor':tumor, 
                'count': count,
                'go':'error'})



    form = FormMutationChoice()       
    return render(request, 'rolls/somatic_interaction_analysis.html', {'form':form,'count': count,})


############ GENE MUTATION ANALYSIS ##############
def gene_mutation_analysis(request):
    count = get_counter()
    if request.method == 'POST':

        form = tumorGeneform(request.POST)
        if form.is_valid(): 
            tumor=request.POST['tumor'] 
            gene=request.POST['gene'] 
            
            inp3=(time.strftime("%Y-%m-%d-%H-%M-%S"))
            dir=os.path.join(output_data, inp3)
            print
            os.makedirs(dir)
            
            out = subprocess.run(['Rscript', 'script/gene_mutation_analysis.R',tumor,gene,dir], capture_output=True, text=True)
            print(out)
            count = increment_counter()
            return_code = out.returncode
            if return_code!=1:
                if os.path.isdir(dir): 
                    files=os.listdir(dir)
                    for file in files:
                        if file[-3:]=='png':
                            if 'lollipopPlot' in file:
                                image_lolli=os.path.join('media/saveanalisi',inp3,file)  
                        if 'txt' in file:
                            result_data=os.path.join(output_data,inp3,file)  
                            result=read_table(result_data)
                    
                    #zip folder analisi -> results.zip
                    folder_to_zip=os.path.join(dir,"results.zip")
                    subprocess.run(["zip", "-r", folder_to_zip, "."],cwd=dir)

                    
                    form=tumorGeneform()
                    return render(request, 'rolls/gene_mutation_analysis.html', {'form':form, 
                        'tumor':tumor,
                        'gene':gene,
                        'dati':result,
                        'image_lolli':image_lolli,
                        'go':'Valid',
                        #'dir':'media/saveanalisi/'+inp3+'/result.txt', #download single table
                        'dir':inp3,
                        'count': count,
                        })

            else:
                    form=tumorGeneform()
                    return render(request, 'rolls/gene_mutation_analysis.html', {'form':form,
                    'gene':gene, 
                    'tumor':tumor,
                    'count': count,
                    'go':'error'})



    form = tumorGeneform()       
    return render(request, 'rolls/gene_mutation_analysis.html', {'form':form,'count': count,})



############ DIFFERENTIAL EXPRESSION MUT VS WT - DESEQ2 #############
def read_table_deseq_demut(file_path):
    txt_data = []

    df = pd.read_csv(file_path, sep='\t',dtype=str) 
    #trasforma padj in valore numerico per sortarlo dal piu piccolo al piu grande
    df['padj'] = pd.to_numeric(df['padj'], errors='coerce')
    df=df[df['padj'] < 0.05].sort_values(by='padj', ascending=True) #.head(500)
    #df = df.sort_values(by='padj', ascending=True).head(500)
    
    df.padj=['%.2E' % Decimal(x) for x in df.padj]
    
   
    # Aggiungi l'intestazione (nomi delle colonne, incluso l'indice) alla lista
    txt_data.append(df.columns.tolist())
    
    # Itera sulle righe del DataFrame e aggiungi ogni riga come lista
    for index, row in df.iterrows():
        txt_data.append(row.tolist())
    return(txt_data)

def de_mut(request):
    count = get_counter()
    if request.method == 'POST':

        form = tumorGeneform(request.POST)
        if form.is_valid(): 
            tumor=request.POST['tumor'] 
            gene=request.POST['gene']
            inp3=(time.strftime("%Y-%m-%d-%H-%M-%S"))
            dir=os.path.join(output_data, inp3)
            
            TCGA_path=get_full_path(config['tcga']['split_count'])

            input_file=os.path.join(TCGA_path,tumor+".tsv")
            os.makedirs(dir)
            
            out = subprocess.run(['Rscript', 'script/MUT_deseq2.R',tumor,gene,dir,input_file], capture_output=True, text=True)
            print(out.returncode)
            count = increment_counter()

            if out.returncode == 1:
                #Errore: non ci sono abbastanza campioni mutati
                form=tumorGeneform()
                return render(request, 'rolls/de_mut.html', {'form':form,
                'tumor':tumor, 
                'gene':gene,
                'count': count,
                'go':'error1'})
                
            elif out.returncode == 2:
                #"Errore: non ci sono abbastanza campioni non mutati.
                form=tumorGeneform()
                return render(request, 'rolls/de_mut.html', {'form':form,
                'tumor':tumor, 
                'gene':gene,
                'count': count,
                'go':'error2'})
                
            elif out.returncode == 0:
                #Lo script R è stato eseguito correttamente
                if os.path.isdir(dir): 
                    files=os.listdir(dir)
                    print(files)
                    n=0
                    for file in files:
                        if 'res' in file:
                            dir_saveresults= os.path.join(output_data, inp3)
                            out_plotly=run([sys.executable,'script/MUT_deseq2_volcano_plotly.py',tumor,file,dir_saveresults],shell=False, stdout=PIPE)
                        
                            #dir=os.path.join('media/saveanalisi',inp3,file)
                            table=os.path.join('media/saveanalisi',inp3,'result.txt')

                            file_txt=os.path.join(output_data,inp3,'result.txt')
                            result=read_table_deseq_demut(file_txt)
                            
                            file_html=tumor+'.html'
                            image_plotly=os.path.join('media/saveanalisi',inp3,file_html)

                        if 'png' in file:
                            
                            if 'Enhanced' in file:
                                image1=os.path.join('media/saveanalisi',inp3,file)  
                                n+=1
                            if 'heatmap' in file:
                                image2=os.path.join('media/saveanalisi',inp3,file)
                                n+=1
                            if 'PCA' in file:
                                image3=os.path.join('media/saveanalisi',inp3,file)
                                n+=1
                            if 'Top50genes' in file:
                                image4=os.path.join('media/saveanalisi',inp3,file)
                                n+=1
                    
                    #zip folder analisi -> results.zip
                    folder_to_zip=os.path.join(dir,"results.zip")
                    subprocess.run(["zip", "-r", folder_to_zip, "."],cwd=dir)
                    if n>1:
                        

                        form=tumorGeneform()
                        return render(request, 'rolls/de_mut.html', {'form':form, 
                            'tumor':tumor,
                            'gene':gene,
                            'image_plotly':image_plotly,
                            'image2':image2,
                            'image3':image3,
                            'image4':image4,
                            'go':'Valid',
                            'dir':inp3,
                            'table':table,
                            'dati':result,
                            'count': count,
                            'parametri':'Mutated gene vs WT'
                            })
                else:
                    form=tumorGeneform()
                    return render(request, 'rolls/de_mut.html', {'form':form,
                    'tumor':tumor, 
                    'gene':gene,
                    'count': count,
                    'go':'error'})

            else:
                form=tumorGeneform()
                return render(request, 'rolls/de_mut.html', {'form':form,
                'tumor':tumor, 
                'gene':gene,
                'count': count,
                'go':'error'})



    form = tumorGeneform()       
    return render(request, 'rolls/de_mut.html', {'form':form,'go':'base','count': count,})


############ DIFFERENTIAL EXPRESSION mutated gene by clinical feature #############
TUMOR_FEATURE_MAPPING_R = {
    "ACC":['gender','person_neoplasm_cancer_status','radiation_therapy'],
    "BLCA":['gender','person_neoplasm_cancer_status','radiation_therapy'],
    "BRCA":['person_neoplasm_cancer_status'], #'radiation_therapy'
    "CESC":['radiation_therapy'],
    "CHOL":['gender','person_neoplasm_cancer_status'],
    "COAD":['gender','person_neoplasm_cancer_status'],#,'radiation_therapy'],
    "DLBC":['gender','radiation_therapy'],
    "ESCA":['alcohol_history_documented','gender','person_neoplasm_cancer_status','radiation_therapy'],
    "GBM":['gender','radiation_therapy'],
    "HNSC":['alcohol_history_documented','gender','person_neoplasm_cancer_status','radiation_therapy'],
    "KICH":['gender','person_neoplasm_cancer_status'],
    "KIRC":['gender','person_neoplasm_cancer_status','radiation_therapy'],
    "KIRP":['gender','person_neoplasm_cancer_status','radiation_therapy'],
    "LGG":['gender','radiation_therapy'],
    "LIHC":['gender','person_neoplasm_cancer_status','radiation_therapy'],
    "LUAD":['gender','person_neoplasm_cancer_status','radiation_therapy'],
    "LUSC":['gender','person_neoplasm_cancer_status','radiation_therapy'],
    "MESO":['gender','person_neoplasm_cancer_status','radiation_therapy'],
    "OV":['person_neoplasm_cancer_status'],
    "PAAD":['alcohol_history_documented','gender','person_neoplasm_cancer_status','radiation_therapy'],
    "PCPG":['gender','radiation_therapy'],
    "PRAD":['radiation_therapy'],
    "READ":['gender','person_neoplasm_cancer_status','radiation_therapy'],
    "SARC":['gender','radiation_therapy'],
    "SKCM":['gender','person_neoplasm_cancer_status','radiation_therapy'],
    "STAD":['gender','person_neoplasm_cancer_status','radiation_therapy'],
    "TGCT":['person_neoplasm_cancer_status','radiation_therapy'],
    "THCA":['gender','person_neoplasm_cancer_status','radiation_therapy'],
    "THYM":['gender','radiation_therapy'],
    "UCEC":['radiation_therapy'],
    "UCS":['radiation_therapy'],
    "UVM":['gender','person_neoplasm_cancer_status','radiation_therapy'],
}

def get_features_R(request):
    # Ottieni il valore del tumore selezionato
    tumor = request.GET.get('tumor')
    
    # Ottieni le feature corrispondenti dal dizionario
    features = TUMOR_FEATURE_MAPPING_R.get(tumor, [])
    
    # Restituisci le feature come tuple (valore, etichetta)
    feature_list = [(feature, feature) for feature in features]
    
    return JsonResponse(feature_list, safe=False)


def de_mut_clinical_feature(request):
    count = get_counter()
    if request.method == 'POST':

        form = featuremutationform(request.POST)
        if form.is_valid(): 
            tumor=request.POST['tumor'] 
            feature=request.POST['feature_selected']
            inp3=(time.strftime("%Y-%m-%d-%H-%M-%S"))
            dir=os.path.join(output_data, inp3)
            
            os.makedirs(dir)
            
            out = subprocess.run(['Rscript', 'script/de_mut_clinical_feature.R',tumor,feature,dir], capture_output=True, text=True)
            print(out)
            count = increment_counter()
            if os.path.isdir(dir): 
                files=os.listdir(dir)
                print(files)
                for file in files:
                    if 'csv' in file:
                        result=file

                    if 'png' in file:
                        if 'forestPlot' in file:
                            image_forest=os.path.join('media/saveanalisi',inp3,file)  
                            
                       
                        if 'coBarplot' in file:
                            image_coBarplot=os.path.join('media/saveanalisi',inp3,file)
                        

                #zip folder analisi -> results.zip
                folder_to_zip=os.path.join(dir,"results.zip")
                subprocess.run(["zip", "-r", folder_to_zip, "."],cwd=dir)


                form=featuremutationform()
                return render(request, 'rolls/de_mut_clinical_feature.html', {'form':form, 
                    'tumor':tumor,
                    'feature':feature,
                    'image_forest':image_forest,
                    'image_coBarplot':image_coBarplot,
                    'go':'Valid',
                    'table':'media/saveanalisi/'+inp3+'/'+result,
                    'dir':inp3,
                    'count': count,
                    })

            else:
                form=featuremutationform()
                return render(request, 'rolls/de_mut_clinical_feature.html', {'form':form,
                'tumor':tumor, 
                'count': count,
                'go':'error'})



    form = featuremutationform()       
    return render(request, 'rolls/de_mut_clinical_feature.html', {'form':form,'count': count,})





############################################
#                                          #
#                CELL TYPES                #
#                                          #
############################################

#analisi:
# 1. Deconvolutio (OK)
# 2. Correlation cell pathway (OK)


########### DECONVOLUTION ########### 
def deconvolution(request):
    count = get_counter()
    if request.method == 'POST':
        form = deconvolution_form(request.POST)
        tumor=request.POST['tumor'] 
        
        dir= os.path.join(base_dir,'deconvolution','results_deconvolution',tumor)
        
        
        if os.path.isdir(dir): 
            inp3=(time.strftime("%Y-%m-%d-%H-%M-%S"))
            dir_saveresults= os.path.join(output_data, inp3)
            print(dir_saveresults)
            os.makedirs(dir_saveresults)

            
            out=run([sys.executable,'script/deconvolution.py',tumor,dir,dir_saveresults],shell=False, stdout=PIPE)
            count = increment_counter()
            # if os.path.isdir(dir): 
            files=os.listdir(dir)
            for file in files:
                print(file)
                if 'boxplot' in file:
                    image_box=os.path.join('media/saveanalisi',inp3,file)
                    print(image_box)
                
                if 'stat_tabel' in file:
                    result_tsv=os.path.join('media/saveanalisi',inp3,file)
                    result_data=os.path.join(output_data,inp3,file)
                    result=read_table(result_data)


            #zip folder analisi -> results.zip
            folder_to_zip=os.path.join(dir_saveresults,"results.zip")
            subprocess.run(["zip", "-r", folder_to_zip, "."],cwd=dir_saveresults)

            form=deconvolution_form()                         
            return render(request, 'rolls/deconvolution.html', {'form':form, 
                'tumor':tumor,
                'image1':image_box,
                'count': count,
                'dati':result,
                'go':'Valid',
                'table':result_tsv,
                'dir':inp3,
                })

        else:
            form=deconvolution_form()
            return render(request, 'rolls/deconvolution.html', {'form':form,
            'count': count,
            'tumor':tumor, 
            'go':'error'})



    form = deconvolution_form()       
    return render(request, 'rolls/deconvolution.html', {'form':form,'count': count,})



########### CORRELATION CELL PATHWAY ########### 
def corr_cell_pathway(request):
    count = get_counter()
    if request.method == 'POST':
        form = formcorrelation(request.POST)
        tumor=request.POST['tumor'] 
        db=request.POST['Db']
        dir= os.path.join(base_dir,'deconvolution','results_correlation_pathways',db,tumor)
        
        
        if os.path.isdir(dir): 
            inp3=(time.strftime("%Y-%m-%d-%H-%M-%S"))
            dir_saveresults= os.path.join(output_data, inp3)
            print(dir_saveresults)
            os.makedirs(dir_saveresults)

            
            out=run([sys.executable,'script/deconvolution.py',tumor,dir,dir_saveresults],shell=False, stdout=PIPE)
            count = increment_counter()
            # if os.path.isdir(dir): 
            files=os.listdir(dir)
            for file in files:
                print(file)
               
                if 'heatmap' in file:
                    image_heat=os.path.join('media/saveanalisi',inp3,file)
                    #print(image_heat)
                if 'tsv' in file:
                    result_tsv=os.path.join('media/saveanalisi',inp3,file)
                    result_data=os.path.join(output_data,inp3,file)
                    result=read_table(result_data)
            
            #zip folder analisi -> results.zip
            folder_to_zip=os.path.join(dir_saveresults,"results.zip")
            subprocess.run(["zip", "-r", folder_to_zip, "."],cwd=dir_saveresults)

            form=formcorrelation()                         
            return render(request, 'rolls/corr_cell_pathway.html', {'form':form, 
                'tumor':tumor,
                'db':db,
                'image2':image_heat,
                'dati':result,
                'go':'Valid',
                #'dir':result_tsv,
                'dir':inp3,
                'count': count,
                })

        else:
            form=formcorrelation()
            return render(request, 'rolls/corr_cell_pathway.html', {'form':form,
            'tumor':tumor,
            'count': count, 
            'go':'error'})



    form = formcorrelation()       
    return render(request, 'rolls/corr_cell_pathway.html', {'form':form,'count': count,})






########################## bozze di analisi non in produzione ##########################

############# CORRELATION ANALYSIS ###################

def correlation_analysis(request):
    count = get_counter()
    if request.method == 'POST':
        form = Analisi_interaction(request.POST)
        if form.is_valid():
            gene=form.cleaned_data['gene']
            miRNA=form.cleaned_data['miRNA']
            tumor= form.cleaned_data['tumor']

            inp3=(time.strftime("%Y-%m-%d-%H-%M-%S"))
            dir= os.path.join(output_data, inp3)
            out=run([sys.executable,'script/overall_survival_interaction.py',miRNA,gene,tumor,dir],shell=False, stdout=PIPE)
            print(out)
            
            dir='rolls/static/media/saveanalisi/'+inp3+'/'
            if os.path.isdir(dir): 
                #dir='rolls/static/media/saveanalisi/'+inp3+'/'
                files=os.listdir(dir)
                images=[]
                for file in files:
                    if file[-3:]=='jpg':
                        images.append('/media/saveanalisi/'+inp3+'/'+file)

                form=Analisi_interaction()
                return render(request, 'rolls/OS_interaction.html', {
                    'form':form, 
                    'formresult': out.stdout.decode('ascii'),
                    'image': images,
                    'go':'Valid',
                    'gene':gene,
                    'tumor':tumor,
                    'miRNA':miRNA,
                    'count': count,
                    })
            else:
                form=Analisi_interaction()
                return render(request, 'rolls/OS_interaction.html', {'form':form, #correlation_analysis
                'miRNA':miRNA,
                'gene':gene,
                'tumor':tumor, 
                'go':'error'})
                

    form=Analisi_interaction()
    return render(request, 'rolls/OS_interaction.html', {'form':form})


def pathwayPROVA(request):
    #menu
    f=open('/mnt/data/notturno/gsva/pathway/name_group_REACTOME.txt').read().rstrip().split('\n')

    #dizionario submenu: submenu2
    tf = open("/mnt/data/notturno/gsva/pathway/Dictionary_name_reactome.json", "r")
    dizionario = json.load(tf)

    #dizionario submenu2: valore msigdb da passare allo script
    dictionary={'chiave1':['v','v2'], 'chiave2':['v1','v2']}


    return render(request, 'rolls/pathwayPROVA.html', {
        'listachiavi':f,
        'dizionary':dictionary,
        'dizionario':dizionario,
    })

#############Overall survival con dati interazione miRNA-mRNA ############

def os_interaction(request):
    if request.method == 'POST':

        if 'interactor' in request.POST:
            form = Gene(request.POST)
            if form.is_valid():
                gene=request.POST['gene']
                tumor=request.POST['tumor']

                out=run([sys.executable,'script/search_interactor.py',gene],shell=False, stdout=PIPE)
                stringa=(out.stdout.decode('ascii')).strip()
                lista=stringa.split(',')
                

                return render(request, 'rolls/OS_interaction.html', {
                    'form':form,
                    'lista':lista,
                    'gene':gene,
                    'tumor':tumor,    
                })
        elif 'Submit' in request.POST:
            go='Selected'

            form = Gene(request.POST)
            gene=request.POST['gene']
            tumor=request.POST['tumor']  

            miRNA=request.POST.get('miRNA', False)

            inp3=(time.strftime("%Y-%m-%d-%H-%M-%S"))

            out=run([sys.executable,'script/overall_survival_interaction.py',gene,miRNA,tumor,inp3],shell=False, stdout=PIPE)
            print(out)
            dir='rolls/static/media/saveanalisi/'+inp3+'/'
            if os.path.isdir(dir): 
                files=os.listdir(dir)
                images=[]
                pvalue=[]
                for file in files:
                    if file[-3:]=='jpg':
                        images.append('/media/saveanalisi/'+inp3+'/'+file)
                pvalue=open(dir+'result.txt').read().rstrip().split("\n")
                
                mylist = zip(images, pvalue)
                
                form=Gene()
            
                return render(request, 'rolls/OS_interaction.html', {
                    'form':form,
                    'gene':gene,
                    'miRNA':miRNA,
                    'tumor':tumor,
                    'go':go,
                    'pvalue':pvalue,
                    'mylist': mylist,
                    'formresult': out.stdout.decode('ascii'),
                    'image': images,
                    'dir':'http://160.80.35.91:7000/static/media/saveanalisi/'+inp3,
                 })
            else:
                form=Gene()
                return render(request, 'rolls/OS_interaction.html', {'form':form,
                'miRNA':miRNA,
                'gene':gene,
                'tumor':tumor, 
                'go':'error'})  

    form = Gene()
    return render(request, 'rolls/OS_interaction.html', {'form':form})
    
