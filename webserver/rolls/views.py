from django.shortcuts import render
from django.http import HttpResponse
from subprocess import run,PIPE
import sys
from matplotlib import image
from rolls.forms import Gene, FormMutationChoice,featuremutationform,formcorrelation,Analisiform_protein,Analisiformcompleto_protein,formSurvival,Analisiform, Deseq2form, Analisiform1, Analisiformcompleto, Analisi_interaction,Analisipath,tumorGeneform,FormTumorMutation
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

from django.http import JsonResponse
from rolls.models import Gene,Pathway,Protein, Gene_symbol

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
    return render(request, 'rolls/home.html')

def documentation(request):
    return render(request, 'rolls/documentation.html')

def read_table(file_path):
    txt_data = []

    df = pd.read_csv(file_path, sep='\t')  # Cambia 'sep' se necessario, es. ',' per CSV
    txt_data.append(df.columns.tolist())
    
    # Itera sulle righe del DataFrame e aggiungi ogni riga come lista
    for index, row in df.iterrows():
        txt_data.append(row.tolist())
    return(txt_data)


def dataset(request):
    file_path = os.path.join(output_data_Table,'table','campioni_TCGA.txt')
    file_path_genetype=os.path.join(output_data_Table,'table','gene_type.txt')
    file_feature_tumor=os.path.join(output_data_Table,'table','Features_tumor.txt')
    file_tumor=os.path.join(output_data_Table,'table','tumor_abbreviations.txt')


    txt_data_clinical=read_table(file_path)

    txt_data_typegene=read_table(file_path_genetype)   
    txt_feature_tumor= read_table(file_feature_tumor)
    txt_tumor= read_table(file_tumor) 

    return render(request, 'rolls/dataset.html', {
        'dati': txt_data_clinical, 
        'dati_genetype':txt_data_typegene,
        'dati_feature':txt_feature_tumor,
        'dati_tumor':txt_tumor,
 })


def contact(request):
    return render(request, 'rolls/contact.html')

def yourdataset(request):
    return render(request,'rolls/yourdataset.html')


##########################################
# Funzioni per l'autocomplete
def gene_suggestions(request):
    if 'term' in request.GET:
        qs = Gene.objects.filter(gene__icontains=request.GET.get('term'))[:10]  # Limita i risultati a 10
        genes = list(qs.values_list('gene', flat=True))
        return JsonResponse(genes, safe=False)
    
def pathway_suggestions(request):
    if 'term' in request.GET:
        qs = Pathway.objects.filter(pathway__icontains=request.GET.get('term'))[:20]  # Limita i risultati a 10
        pathways = list(qs.values_list('pathway', flat=True))
        return JsonResponse(pathways, safe=False)
    
def protein_suggestions(request):
    if 'term' in request.GET:
        qs = Protein.objects.filter(protein__icontains=request.GET.get('term'))[:20]  # Limita i risultati a 10
        protein = list(qs.values_list('protein', flat=True))
        return JsonResponse(protein, safe=False)
    
def gene_symbol_suggestions(request):
    if 'term' in request.GET:
        qs = Gene_symbol.objects.filter(gene_symbol__icontains=request.GET.get('term'))[:20]  # Limita i risultati a 10
        gene_symbol = list(qs.values_list('gene_symbol', flat=True))
        return JsonResponse(gene_symbol, safe=False)
##########################################




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







######### OVERALL SURVIVAL ################
def overall_survival(request):
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
            

            if os.path.isdir(dir): 
                files=os.listdir(dir)
                for file in files:
                    if 'jpeg' in file:
                        image=os.path.join('media/saveanalisi',inp3,file)
                        
                        form=formSurvival()

                        return render(request, 'rolls/overall_survival.html', {
                            'form':form, 
                            'formresult': out.stdout.decode('ascii'),
                            'image': image,
                            'go':'Valid',
                            'gene':gene ,
                            'tumor':tumor,
                            'method':methods,
                            'dir':inp3})

            else:
                form=formSurvival()
                return render(request, 'rolls/overall_survival.html', {'form':form,
                'formresult': out.stdout.decode('ascii'),
                'gene':gene,
                'tumor':tumor, 
                'go':'error'})


    form=formSurvival()
    return render(request, 'rolls/overall_survival.html', {'form':form})


###########Overall survival con pathway activity score######
def os_pathway(request):
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
            
            dir='rolls/static/media/saveanalisi/'+inp3+'/'
            if os.path.isdir(dir): 
                files=os.listdir(dir)
                for file in files:
                    if file[-3:]=='png':
                        image='/media/saveanalisi/'+inp3+'/'+file
                form=Analisipath()
                return render(request, 'rolls/OS_pathway.html', {
                    'form':form, 
                    'formresult': out.stdout.decode('ascii'),
                    'image': image,
                    'go':'Valid',
                    'pathway':pathway ,
                    'tumor':tumor,
                    'method':method,
                    })
            else:
                form=Analisipath()
                return render(request, 'rolls/OS_pathway.html', {'form':form,
                'pathway':pathway ,
                'tumor':tumor, 
                'go':'error'})


    form=Analisipath()
    return render(request, 'rolls/OS_pathway.html', {'form':form})


######### survival with gene mutation status ######### 
def survival_with_gene_mutation_status(request):
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
            
            if os.path.isdir(dir): 
                files=os.listdir(dir)
                for file in files:
                    if 'png' in file:
                        image=os.path.join('media/saveanalisi',inp3,file)    
                
                form=tumorGeneform()
                return render(request, 'rolls/survival_with_gene_mutation_status.html', {'form':form, 
                    'gene':gene,
                    'tumor':tumor,
                    'image':image,
                    'go':'Valid',
                    'dir':inp3,
                    })

            else:
                
                form=tumorGeneform()
                return render(request, 'rolls/survival_with_gene_mutation_status.html', {
                'form':form,
                'formresult':'analysis is not available for the entered name',                                                                      
                'tumor':tumor, 
                'go':'error'})



    form = tumorGeneform()       
    return render(request, 'rolls/survival_with_gene_mutation_status.html', {'form':form})




############     DIFFERENTIAL EXPRESSION SINGLE TUMOR  TRASCRITTOMIC   ############ manca da gestire errori
def diff_exp_single_tumor(request):
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
                print(out)
                
                
                if os.path.isdir(dir): 
                    files=os.listdir(dir)
                    n=0
                    for file in files:
                        
                        if file[-4:]=='html':
                            n+=1
                           # image='/media/saveanalisi/'+time_dir+'/'+file
                            image=os.path.join('media/saveanalisi',time_dir,file)
                            form=Analisiformcompleto()
                            return render(request, 'rolls/diff_exp_single_tumor.html', {
                                'form':form, 
                                'formresult': out.stdout.decode('ascii'),
                                'image': image,
                                'go':'Valid',
                                'gene':gene,
                                'tumor':tumor,
                                'feature':feature,
                                'parametri': parametri[feature],})
                    if n==0:
                      
                        form=Analisiformcompleto()
                        return render(request, 'rolls/diff_exp_single_tumor.html', {'form':form,
                        'formresult': out,
                        'feature': feature,
                        'tumor':tumor, 
                        'go':'error'})
                else:
                    form=Analisiformcompleto()
                    return render(request, 'rolls/diff_exp_single_tumor.html', {'form':form,
                    'formresult': out.stdout.decode('ascii'),
                    'feature': feature,
                    'tumor':tumor, 
                    'go':'error'})


    form=Analisiformcompleto()
    return render(request, 'rolls/diff_exp_single_tumor.html', {'form':form})



#############      DIFFERENTIAL EXPRESSION ANALYSIS ALL TUMOR FOR FEATURE TRASCRITTOMIC   ############# manca da gestire errori
def differential_expression(request):
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
            
            if os.path.isdir(dir): 
                files=os.listdir(dir)
                for file in files:
                    if file[-3:]=='jpg':
                        image=os.path.join('media/saveanalisi',inp3,file)
                        #image='media/saveanalisi/'+inp3+'/'+file

                    if 'result' in file:
                        result_data=os.path.join(output_data,inp3,file)
                        result=read_table(result_data)
                    
                form=Analisiform()
                return render(request, 'rolls/differential_expression.html', {
                    'form':form, 
                    'formresult': out.stdout.decode('ascii'),
                    'image': image,
                    'go':'Valid',
                    'gene':gene,
                    'feature':feature,
                    'parametri': parametri[feature],
                    'dir':'media/saveanalisi/'+inp3+'/result.txt',
                    'dati':result,
                    })
            else:
                form=Analisiform()
                return render(request, 'rolls/differential_expression.html', {
                    'form':form,
                    'formresult': out.stdout.decode('ascii'),
                    'feature': feature,
                    'gene':gene, 
                    'go':'error'})


    form=Analisiform()
    return render(request, 'rolls/differential_expression.html', {'form':form})






#############      DIFFERENTIAL EXPRESSION ANALYSIS ALL TUMOR FOR FEATURE   PROTEOMIC   ############# manca da gestire errori
def differential_expression_protein(request):
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
            
            if os.path.isdir(dir): 
                files=os.listdir(dir)
                for file in files:
                    if file[-3:]=='jpg':
                        image=os.path.join('media/saveanalisi',inp3,file)
                        #image='media/saveanalisi/'+inp3+'/'+file
                    if 'txt' in file:
                        result_data=os.path.join(output_data,inp3,file)
                        result=read_table(result_data)
                form=Analisiform_protein()
                return render(request, 'rolls/differential_expression_protein.html', {
                    'form':form, 
                    'formresult': out.stdout.decode('ascii'),
                    'image': image,
                    'go':'Valid',
                    'gene':gene,
                    'feature':feature,
                    'parametri': parametri[feature],
                    'dir':'media/saveanalisi/'+inp3+'/result.txt',
                    'dati':result,
                    })
            else:
                form=Analisiform_protein()
                return render(request, 'rolls/differential_expression_protein.html', {
                    'form':form,
                    'formresult': out.stdout.decode('ascii'),
                    'feature': feature,
                    'gene':gene, 
                    'go':'error'})


    form=Analisiform_protein()
    return render(request, 'rolls/differential_expression_protein.html', {'form':form})


############     DIFFERENTIAL EXPRESSION SINGLE TUMOR PROTEOMIC   ############ manca da gestire errori
def diff_exp_single_tumor_protein(request):
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
                
                
                if os.path.isdir(dir): 
                    files=os.listdir(dir)
                    for file in files:
                        if file[-4:]=='html':
                           # image='/media/saveanalisi/'+time_dir+'/'+file
                            image=os.path.join('media/saveanalisi',time_dir,file)
                            form=Analisiformcompleto_protein()
                            return render(request, 'rolls/diff_exp_single_tumor_protein.html', {
                                'form':form, 
                                'formresult': out.stdout.decode('ascii'),
                                'image': image,
                                'go':'Valid',
                                'gene':gene,
                                'tumor':tumor,
                                'feature':feature,
                                'parametri': parametri[feature],})
                else:
                    form=Analisiformcompleto_protein()
                    return render(request, 'rolls/diff_exp_single_tumor_protein.html', {'form':form,
                    'feature': feature,
                    'tumor':tumor, 
                    'go':'error'})


    form=Analisiformcompleto_protein()
    return render(request, 'rolls/diff_exp_single_tumor_protein.html', {'form':form})




####### DESEQ2 analisi###############

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
    

    df = pd.read_csv(file_path, sep='\t')  # Cambia 'sep' se necessario, es. ',' per CSV
    df = df.sort_values(by='padj', ascending=True).head(500)

    
   
    # Aggiungi l'intestazione (nomi delle colonne, incluso l'indice) alla lista
    txt_data.append(df.columns.tolist())
    
    # Itera sulle righe del DataFrame e aggiungi ogni riga come lista
    for index, row in df.iterrows():
        txt_data.append(row.tolist())
    return(txt_data)


def deseq2(request):
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
                    'dir':'media/saveanalisi/'+inp3+'/'+result_file, 
                    'dati':result,
                    })

            else:
                form=Deseq2form()
                return render(request, 'rolls/deseq2.html', {'form':form,
                'feature': feature,
                'tumor':tumor, 
                'go':'error'})



    form = Deseq2form()       
    return render(request, 'rolls/deseq2.html', {'form':form})


# def deseq2_old(request):
#     if request.method == 'POST':

#         if 'features' in request.POST:
#             form = Deseq2form(request.POST)
#             if form.is_valid():
#                 tumor=request.POST['tumor'] 
#                 out=run([sys.executable,'script/search_feature_deseq2.py',tumor],shell=False, stdout=PIPE)
#                 stringa=(out.stdout.decode('ascii')).strip()
#                 lista=stringa.split(',')

#                 return render(request, 'rolls/deseq2.html', {
#                     'form':form,
#                     'lista':lista,
#                     'tumor':tumor,
#                     })
#         elif 'Submit' in request.POST:
#             go='Valid'
#             form = Deseq2form(request.POST)
#             tumor=request.POST['tumor'] 
#             feature=request.POST.get('feature', False)
#             dir= os.path.join(base_dir,'deseq2', feature,tumor)
            
            
#             if os.path.isdir(dir): 
#                 inp3=(time.strftime("%Y-%m-%d-%H-%M-%S"))
#                 dir_saveresults= os.path.join(output_data, inp3)
#                 print(dir_saveresults)
#                 os.makedirs(dir_saveresults)

#                 result_file='result_' + tumor + '_2.txt'
#                 out=run([sys.executable,'script/deseq2.py',tumor,dir,dir_saveresults,result_file],shell=False, stdout=PIPE)
                
#                 images=choseimage(tumor,dir,dir_saveresults)
                
#                 file_txt=os.path.join(output_data,inp3,result_file)
#                 print(file_txt)
#                 result=read_table_deseq(file_txt)
#                 form=Deseq2form()                         
#                 return render(request, 'rolls/deseq2.html', {'form':form, 
#                     'feature': feature,
#                     'tumor':tumor,
#                     'enhancedimage': os.path.join('media/saveanalisi',inp3,images[0]),
#                     'images1': os.path.join('media/saveanalisi',inp3,images[1]),
#                     'images2': os.path.join('media/saveanalisi',inp3,images[2]),
#                     'images3': os.path.join('media/saveanalisi',inp3,images[3]),
#                     'image_plotly':os.path.join('media/saveanalisi',inp3,images[4]),
#                     'go':'Valid',
#                     'parametri': parametri[feature],
#                     'dir':'media/saveanalisi/'+inp3+'/'+result_file, 
#                     'dati':result,
#                     })

#             else:
#                 form=Deseq2form()
#                 return render(request, 'rolls/deseq2.html', {'form':form,
#                 'feature': feature,
#                 'tumor':tumor, 
#                 'go':'error'})



#     form = Deseq2form()       
#     return render(request, 'rolls/deseq2.html', {'form':form})




########### MUTATION ANALYSES ########### 
def tumor_mutation_analysis(request):
    if request.method == 'POST':

        form = FormTumorMutation(request.POST)
        if form.is_valid(): 
            tumor=request.POST['tumor'] 
            
            inp3=(time.strftime("%Y-%m-%d-%H-%M-%S"))
            dir=os.path.join(output_data, inp3)
            
            os.makedirs(dir)
            
            out = subprocess.run(['Rscript', 'script/tumor_mutation_analysis.R',tumor,dir], capture_output=True, text=True)
            print(out)
            
            if os.path.isdir(dir): 
                files=os.listdir(dir)
                print(files)
                for file in files:
                    if file[-3:]=='png':
                        if 'summary' in file:
                            image_summary=os.path.join('media/saveanalisi',inp3,file)  
                            
                       
                        if 'Titv' in file:
                            image_titv=os.path.join('media/saveanalisi',inp3,file)
                        

                        
                form=FormTumorMutation()
                return render(request, 'rolls/tumor_mutation_analysis.html', {'form':form, 
                    'tumor':tumor,
                    'image_summary':image_summary,
                    'image_titv':image_titv,
                    'go':'Valid',
                    'dir':inp3,
                    })

            else:
                form=FormTumorMutation()
                return render(request, 'rolls/tumor_mutation_analysis.html', {'form':form,
                'tumor':tumor, 
                'go':'error'})



    form = FormTumorMutation()       
    return render(request, 'rolls/tumor_mutation_analysis.html', {'form':form})


############## oncoplot #############
def tumor_oncoplot(request):
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
            
            if os.path.isdir(dir): 
                files=os.listdir(dir)
                print(files)
                for file in files:
                    if file[-3:]=='png':                            
                        if 'oncoplot' in file:
                            image_oncoplot=os.path.join('media/saveanalisi',inp3,file)
                    
                        
                form=FormMutationChoice()
                return render(request, 'rolls/tumor_oncoplot.html', {'form':form, 
                    'tumor':tumor,
                    'image_oncoplot':image_oncoplot,
                    'go':'Valid',
                    'dir':inp3,
                    })

            else:
                form=FormMutationChoice()
                return render(request, 'rolls/tumor_oncoplot.html', {'form':form,
                'tumor':tumor, 
                'go':'error'})



    form = FormMutationChoice()       
    return render(request, 'rolls/tumor_oncoplot.html', {'form':form})



###########differential expression mutato vs non mutato ############# da implementare
def de_mut(request):
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
            print(out)
            
            if os.path.isdir(dir): 
                files=os.listdir(dir)
                print(files)
                n=0
                for file in files:
                    if 'res' in file:
                        dir=os.path.join('media/saveanalisi',inp3,file)
                        file_txt=os.path.join(output_data,inp3,file)
                        result=read_table_deseq(file_txt)
                        #dir_saveresults= os.path.join(output_data, inp3)
                        # out_plotly=run([sys.executable,'script/MUT_deseq2.py',tumor,file_txt,dir_saveresults],shell=False, stdout=PIPE)

                        # if 'html' in file:
                        #     image_plotly=file


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

                    
                if n>1:
                    form=tumorGeneform()
                    return render(request, 'rolls/de_mut.html', {'form':form, 
                        'tumor':tumor,
                        'gene':gene,
                        #'image_plotly':image_plotly,
                        'image1':image1,
                        'image2':image2,
                        'image3':image3,
                        'image4':image4,
                        'go':'Valid',
                        'dir':dir,
                        'dati':result,
                        })
                else:
                    form=tumorGeneform()
                    return render(request, 'rolls/de_mut.html', {'form':form,
                    'tumor':tumor, 
                    'gene':gene,
                    'go':'error'})

            else:
                form=tumorGeneform()
                return render(request, 'rolls/de_mut.html', {'form':form,
                'tumor':tumor, 
                'gene':gene,
                'go':'error'})



    form = tumorGeneform()       
    return render(request, 'rolls/de_mut.html', {'form':form,'go':'base'})


#DE_mutated gene by clinical feature
TUMOR_FEATURE_MAPPING_R = {
    "ACC":['gender','person_neoplasm_cancer_status','radiation_therapy'],
    "BLCA":['gender','person_neoplasm_cancer_status','radiation_therapy'],
    "BRCA":['person_neoplasm_cancer_status','radiation_therapy'],
    "CESC":['radiation_therapy'],
    "CHOL":['gender','person_neoplasm_cancer_status'],
    "COAD":['gender','person_neoplasm_cancer_status','radiation_therapy'],
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
                        

                        
                form=featuremutationform()
                return render(request, 'rolls/de_mut_clinical_feature.html', {'form':form, 
                    'tumor':tumor,
                    'feature':feature,
                    'image_forest':image_forest,
                    'image_coBarplot':image_coBarplot,
                    'go':'Valid',
                    'dir':'media/saveanalisi/'+inp3+'/'+result,
                    })

            else:
                form=featuremutationform()
                return render(request, 'rolls/de_mut_clinical_feature.html', {'form':form,
                'tumor':tumor, 
                'go':'error'})



    form = featuremutationform()       
    return render(request, 'rolls/de_mut_clinical_feature.html', {'form':form})



############somatic_interaction_analysis ########

def read_table_comma(file_path):
    txt_data = []

    df = pd.read_csv(file_path, sep=',')  # Cambia 'sep' se necessario, es. ',' per CSV
    txt_data.append(df.columns.tolist())
    
    # Itera sulle righe del DataFrame e aggiungi ogni riga come lista
    for index, row in df.iterrows():
        txt_data.append(row.tolist())
    return(txt_data)

def somatic_interaction_analysis(request):
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
                        
                form=FormMutationChoice()
                return render(request, 'rolls/somatic_interaction_analysis.html', {'form':form, 
                    'tumor':tumor,
                    'number':number,
                    'image_interact':image_interact,
                    'go':'Valid',
                    'dir':inp3,
                    'dati':result,
                    'dir':'media/saveanalisi/'+inp3+'/'+name,
                    })

            else:
                form=FormMutationChoice()
                return render(request, 'rolls/somatic_interaction_analysis.html', {'form':form,
                'tumor':tumor, 
                'go':'error'})



    form = FormMutationChoice()       
    return render(request, 'rolls/somatic_interaction_analysis.html', {'form':form})



############ gene mutation analysis ##############
def gene_mutation_analysis(request):
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
            return_code = out.returncode
            if return_code!=1:
                if os.path.isdir(dir): 
                    files=os.listdir(dir)
                    for file in files:
                        if file[-3:]=='png':
                            # if 'somaticInteractions' in file:
                            #     image_interact=os.path.join('media/saveanalisi',inp3,file)  
                            # if 'TumorVAF' in file: 
                            #     image_VAF=os.path.join('media/saveanalisi',inp3,file)  
                            if 'lollipopPlot' in file:
                                image_lolli=os.path.join('media/saveanalisi',inp3,file)  
                        if 'txt' in file:
                            result_data=os.path.join(output_data,inp3,file)  
                            result=read_table(result_data)
                    form=tumorGeneform()
                    return render(request, 'rolls/gene_mutation_analysis.html', {'form':form, 
                        'tumor':tumor,
                        'gene':gene,
                        'dati':result,
                        'image_lolli':image_lolli,
                        'go':'Valid',
                        'dir':'media/saveanalisi/'+inp3+'/result.txt',
                        })

            else:
                    form=tumorGeneform()
                    return render(request, 'rolls/gene_mutation_analysis.html', {'form':form,
                    'gene':gene, 
                    'tumor':tumor,
                    'go':'error'})



    form = tumorGeneform()       
    return render(request, 'rolls/gene_mutation_analysis.html', {'form':form})



########### DECONVOLUTION ########### 
def deconvolution(request):
    if request.method == 'POST':
        form = Deseq2form(request.POST)
        tumor=request.POST['tumor'] 
        
        dir= os.path.join(base_dir,'deconvolution','results_deconvolution',tumor)
        
        
        if os.path.isdir(dir): 
            inp3=(time.strftime("%Y-%m-%d-%H-%M-%S"))
            dir_saveresults= os.path.join(output_data, inp3)
            print(dir_saveresults)
            os.makedirs(dir_saveresults)

            
            out=run([sys.executable,'script/deconvolution.py',tumor,dir,dir_saveresults],shell=False, stdout=PIPE)
            
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
                    
            form=Deseq2form()                         
            return render(request, 'rolls/deconvolution.html', {'form':form, 
                'tumor':tumor,
                'image1':image_box,
                
                'dati':result,
                'go':'Valid',
                'dir':result_tsv,
                })

        else:
            form=Deseq2form()
            return render(request, 'rolls/deconvolution.html', {'form':form,
            'tumor':tumor, 
            'go':'error'})



    form = Deseq2form()       
    return render(request, 'rolls/deconvolution.html', {'form':form})


def corr_cell_pathway(request):
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
                    
            form=formcorrelation()                         
            return render(request, 'rolls/corr_cell_pathway.html', {'form':form, 
                'tumor':tumor,
                'db':db,
                'image2':image_heat,
                'dati':result,
                'go':'Valid',
                'dir':result_tsv,
                })

        else:
            form=formcorrelation()
            return render(request, 'rolls/corr_cell_pathway.html', {'form':form,
            'tumor':tumor, 
            'go':'error'})



    form = formcorrelation()       
    return render(request, 'rolls/corr_cell_pathway.html', {'form':form})


#### bozze da revisionare##########################

############# CORRELATION ANALYSIS ###################

def correlation_analysis(request):
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
    
