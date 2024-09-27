from django.shortcuts import render
from django.http import HttpResponse
from subprocess import run,PIPE
import sys
from matplotlib import image
from rolls.forms import Gene, Analisiform_protein,Analisiformcompleto_protein,formSurvival,Analisiform, Deseq2form, Analisiform1, Analisiformcompleto, Analisi_interaction,Analisipath,tumorGeneform,FormTumorMutation
import os
from django.http import StreamingHttpResponse
from wsgiref.util import FileWrapper
import mimetypes
from shutil import make_archive
import time
import os.path
from django.conf import settings
import shutil

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

def dataset(request):
    return render(request, 'rolls/dataset.html')


def contact(request):
    return render(request, 'rolls/contact.html')

def yourdataset(request):
    return render(request,'rolls/yourdataset.html')

####### autocomplete ##############
#funzionante prova1:
# def gene_suggestions(request):
#     if 'term' in request.GET:
#         query = request.GET.get('term')
#         genes = Gene.objects.filter(gene__icontains=query)[:10]  # Limita a 10 risultati
#         gene_names = list(genes.values_list('gene', flat=True))
#         return JsonResponse(gene_names, safe=False)
#     return JsonResponse([], safe=False)


# Funzione per l'autocomplete prova2
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
            os.makedirs(dir)
            out=run([sys.executable,'script/overall_survival.py',gene,tumor,dir,methods],shell=False, stdout=PIPE)
            print(out)
            

            if os.path.isdir(dir): 
                files=os.listdir(dir)
                for file in files:
                    if file[-3:]=='png':
                        image=os.path.join('media/saveanalisi',inp3,file)
                        
                form=formSurvival()

                return render(request, 'rolls/overall_survival.html', {
                    'form':form, 
                    #'formresult': out.stdout.decode('ascii'),
                    'image': image,
                    'go':'Valid',
                    'gene':gene ,
                    'tumor':tumor,
                    'method':methods,
                    'dir':inp3})
            else:
                form=formSurvival()
                return render(request, 'rolls/overall_survival.html', {'form':form,
                'gene':gene,
                'tumor':tumor, 
                'go':'error'})


    form=formSurvival()
    return render(request, 'rolls/overall_survival.html', {'form':form})


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


############     DIFFERENTIAL EXPRESSION SINGLE TUMOR  TRASCRITTOMIC   ############
def diff_exp_single_tumor(request):
    if request.method == 'POST':
            form = Analisiformcompleto(request.POST)
            if form.is_valid():
                gene=form.cleaned_data['gene']
                tumor=form.cleaned_data['tumor']
                feature=form.cleaned_data['feature']
                
                time_dir=(time.strftime("%Y-%m-%d-%H-%M-%S"))
                dir= os.path.join(output_data, time_dir)
               # dir=os.path.join(settings.BASE_DIR, 'rolls', 'static', 'media', 'saveanalisi', time_dir)
                os.makedirs(dir)
                out=run([sys.executable,'script/Differential_expression_boxplot_plotly_new.py',gene,tumor,feature,dir],shell=False, stdout=PIPE)
                print(out)
                
                
                if os.path.isdir(dir): 
                    files=os.listdir(dir)
                    for file in files:
                        if file[-4:]=='html':
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
                else:
                    form=Analisiformcompleto()
                    return render(request, 'rolls/diff_exp_single_tumor.html', {'form':form,
                    'feature': feature,
                    'tumor':tumor, 
                    'go':'error'})


    form=Analisiformcompleto()
    return render(request, 'rolls/diff_exp_single_tumor.html', {'form':form})

#############      DIFFERENTIAL EXPRESSION ANALYSIS ALL TUMOR FOR FEATURE TRASCRITTOMIC   #############
def differential_expression(request):
    if request.method == 'POST':
        form = Analisiform(request.POST)
        if form.is_valid():
            gene=form.cleaned_data['gene']
            feature=form.cleaned_data['feature']
            
            inp3=(time.strftime("%Y-%m-%d-%H-%M-%S"))
            
            dir= os.path.join(output_data, inp3)
            out=run([sys.executable,'script/boxplot_all_tumor_giusto_new.py',gene,feature,dir],shell=False, stdout=PIPE)
            print(out)
            
            if os.path.isdir(dir): 
                files=os.listdir(dir)
                for file in files:
                    if file[-3:]=='jpg':
                        image=os.path.join('media/saveanalisi',inp3,file)
                        #image='media/saveanalisi/'+inp3+'/'+file
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





#############      DIFFERENTIAL EXPRESSION ANALYSIS ALL TUMOR FOR FEATURE   PROTEOMIC   #############
def differential_expression_protein(request):
    if request.method == 'POST':
        form = Analisiform_protein(request.POST)
        if form.is_valid():
            gene=form.cleaned_data['gene']
            feature=form.cleaned_data['feature']
            
            inp3=(time.strftime("%Y-%m-%d-%H-%M-%S"))
            
            dir= os.path.join(output_data, inp3)
            out=run([sys.executable,'script/boxplot_all_tumor_giusto_new.py',gene,feature,dir],shell=False, stdout=PIPE)
            print(out)
            
            if os.path.isdir(dir): 
                files=os.listdir(dir)
                for file in files:
                    if file[-3:]=='jpg':
                        image=os.path.join('media/saveanalisi',inp3,file)
                        #image='media/saveanalisi/'+inp3+'/'+file
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


############     DIFFERENTIAL EXPRESSION SINGLE TUMOR PROTEOMIC   ############
def diff_exp_single_tumor_protein(request):
    if request.method == 'POST':
            form = Analisiformcompleto_protein(request.POST)
            if form.is_valid():
                gene=form.cleaned_data['gene']
                tumor=form.cleaned_data['tumor']
                feature=form.cleaned_data['feature']
                protein='protein'
                time_dir=(time.strftime("%Y-%m-%d-%H-%M-%S"))
                dir= os.path.join(output_data, time_dir)
               # dir=os.path.join(settings.BASE_DIR, 'rolls', 'static', 'media', 'saveanalisi', time_dir)
                os.makedirs(dir)
                out=run([sys.executable,'script/Differential_expression_boxplot_plotly_new.py',gene,tumor,feature,dir,protein],shell=False, stdout=PIPE)
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


def deseq2(request):
    if request.method == 'POST':

        if 'features' in request.POST:
            form = Deseq2form(request.POST)
            if form.is_valid():
                tumor=request.POST['tumor'] 
                out=run([sys.executable,'script/search_feature_deseq2.py',tumor],shell=False, stdout=PIPE)
                stringa=(out.stdout.decode('ascii')).strip()
                lista=stringa.split(',')

                return render(request, 'rolls/deseq2.html', {
                    'form':form,
                    'lista':lista,
                    'tumor':tumor,
                    })
        elif 'Submit' in request.POST:
            go='Valid'
            form = Deseq2form(request.POST)
            tumor=request.POST['tumor'] 
            feature=request.POST.get('feature', False)
            dir= os.path.join(base_dir,'deseq2', feature,tumor)
            
            
            if os.path.isdir(dir): 
                inp3=(time.strftime("%Y-%m-%d-%H-%M-%S"))
                dir_saveresults= os.path.join(output_data, inp3)
                print(dir_saveresults)
                os.makedirs(dir_saveresults)

                
                out=run([sys.executable,'script/deseq2.py',tumor,dir,dir_saveresults],shell=False, stdout=PIPE)
                
                images=choseimage(tumor,dir,dir_saveresults)
                
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
                    'dir':'media/saveanalisi/'+inp3+'/result_' + tumor + '.txt' 
                    })

            else:
                form=Deseq2form()
                return render(request, 'rolls/deseq2.html', {'form':form,
                'feature': feature,
                'tumor':tumor, 
                'go':'error'})



    form = Deseq2form()       
    return render(request, 'rolls/deseq2.html', {'form':form})



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
                            
                        if 'oncoplot' in file:
                            image_oncoplot=os.path.join('media/saveanalisi',inp3,file)
                        if 'Titv' in file:
                            image_titv=os.path.join('media/saveanalisi',inp3,file)

                        
                form=FormTumorMutation()
                return render(request, 'rolls/tumor_mutation_analysis.html', {'form':form, 
                    'tumor':tumor,
                    'image_summary':image_summary,
                    'image_oncoplot':image_oncoplot,
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
            
            if os.path.isdir(dir): 
                files=os.listdir(dir)
                for file in files:
                    if file[-3:]=='png':
                        if 'somaticInteractions' in file:
                            image_interact=os.path.join('media/saveanalisi',inp3,file)  
                        if 'TumorVAF' in file: 
                            image_VAF=os.path.join('media/saveanalisi',inp3,file)  
                        if 'lollipopPlot' in file:
                            image_lolli=os.path.join('media/saveanalisi',inp3,file)  

                form=tumorGeneform()
                return render(request, 'rolls/gene_mutation_analysis.html', {'form':form, 
                    'tumor':tumor,
                    'gene':gene,
                    'image_interact':image_interact,
                    'image_VAF':image_VAF,
                    'image_lolli':image_lolli,
                    'go':'Valid',
                    'dir':inp3,
                    })

            else:
                    form=tumorGeneform()
                    return render(request, 'rolls/gene_mutation_analysis.html', {'form':form,
                    'gene':gene, 
                    'tumor':tumor,
                    'go':'error'})



    form = tumorGeneform()       
    return render(request, 'rolls/gene_mutation_analysis.html', {'form':form})




def survival_with_gene_mutation_status(request):
    if request.method == 'POST':

        form = Deseq2form(request.POST)
        if form.is_valid(): 
            tumor=request.POST['tumor'] 
            
            inp3=(time.strftime("%Y-%m-%d-%H-%M-%S"))
            dir=os.path.join(output_data, inp3)
            print
            os.makedirs(dir)
            
            out = subprocess.run(['Rscript', 'script/survival_with_gene_mutation_status.R',tumor,dir], capture_output=True, text=True)
            print(out)
            
            if os.path.isdir(dir): 
                files=os.listdir(dir)
                for file in files:
                    if file[-3:]=='png':
                        image=os.path.join('media/saveanalisi',inp3,file)    
                print(image)
                form=Deseq2form()
                return render(request, 'rolls/survival_with_gene_mutation_status.html', {'form':form, 
            
                    'tumor':tumor,
                    'image':image,
                    'go':'Valid',
                    'dir':inp3,
                    })

            else:
                    form=Deseq2form()
                    return render(request, 'rolls/survival_with_gene_mutation_status.html', {'form':form,
                    'tumor':tumor, 
                    'go':'error'})



    form = Deseq2form()       
    return render(request, 'rolls/survival_with_gene_mutation_status.html', {'form':form})



########### DECONVOLUTION ########### 
def deconvolution(request):
    if request.method == 'POST':
        form = Deseq2form(request.POST)
        tumor=request.POST['tumor'] 
        
        dir= os.path.join(base_dir,'deconvolution','results',tumor)
        
        
        if os.path.isdir(dir): 
            inp3=(time.strftime("%Y-%m-%d-%H-%M-%S"))
            dir_saveresults= os.path.join(output_data, inp3)
            print(dir_saveresults)
            os.makedirs(dir_saveresults)

            
            out=run([sys.executable,'script/deconvolution.py',tumor,dir,dir_saveresults],shell=False, stdout=PIPE)
            
            if os.path.isdir(dir): 
                files=os.listdir(dir)
                for file in files:
                    print(file)
                    if 'jpeg' in file:
                        image=os.path.join('media/saveanalisi',inp3,file)
                        print(image)
                    if 'tsv' in file:
                        result_tsv=os.path.join('media/saveanalisi',inp3,file)

                    
            form=Deseq2form()                         
            return render(request, 'rolls/deconvolution.html', {'form':form, 
                'tumor':tumor,
                'image':image,
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



