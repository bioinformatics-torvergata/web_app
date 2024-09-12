from django.shortcuts import render
from django.http import HttpResponse
from subprocess import run,PIPE
import sys
from matplotlib import image
from rolls.forms import Gene, Analisiform, Deseq2form, Analisiform1, Analisiformcompleto, Analisi_interaction,Analisipath
import os
from django.http import StreamingHttpResponse
from wsgiref.util import FileWrapper
import mimetypes
from shutil import make_archive
import time
import os.path
from django.conf import settings


from django.http import JsonResponse
from rolls.models import Gene
from rolls.models import Pathway


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
        form = Analisiform1(request.POST)
        if form.is_valid():
            gene=form.cleaned_data['gene']
            tumor=form.cleaned_data['tumor']
            print(gene, tumor)
            
            inp3=(time.strftime("%Y-%m-%d-%H-%M-%S"))
            out=run([sys.executable,'script/overall_survival.py',gene,tumor,inp3],shell=False, stdout=PIPE)
            print(out)
            dir= os.path.join(output_data, inp3)
            #dir='rolls/static/media/saveanalisi/'+inp3+'/'
            if os.path.isdir(dir): 
                files=os.listdir(dir)
                for file in files:
                    if file[-3:]=='png':
                        image=os.path.join('media/saveanalisi',inp3,file)
                        #image='/media/saveanalisi/'+inp3+'/'+file
                form=Analisiform1()

                return render(request, 'rolls/overall_survival.html', {
                    'form':form, 
                    'formresult': out.stdout.decode('ascii'),
                    'image': image,
                    'go':'Valid',
                    'gene':gene ,
                    'tumor':tumor,
                    'dir':inp3})
            else:
                form=Analisiform1()
                return render(request, 'rolls/overall_survival.html', {'form':form,
                'gene':gene,
                'tumor':tumor, 
                'go':'error'})


    form=Analisiform1()
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
            
            inp3=(time.strftime("%Y-%m-%d-%H-%M-%S"))

            out=run([sys.executable,'script/OS_pathway.py',tumor,pathway,inp3],shell=False, stdout=PIPE)
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
                    'tumor':tumor,})
            else:
                form=Analisipath()
                return render(request, 'rolls/OS_pathway.html', {'form':form,
                'pathway':pathway ,
                'tumor':tumor, 
                'go':'error'})


    form=Analisipath()
    return render(request, 'rolls/OS_pathway.html', {'form':form})


############     DIFFERENTIAL EXPRESSION SINGLE TUMOR     ############
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



#############      DIFFERENTIAL EXPRESSION ANALYSIS ALL TUMOR FOR FEATURE      #############
def differential_expression(request):
    if request.method == 'POST':
        form = Analisiform(request.POST)
        if form.is_valid():
            gene=form.cleaned_data['gene']
            feature=form.cleaned_data['feature']
            
            inp3=(time.strftime("%Y-%m-%d-%H-%M-%S"))
            
            dir= os.path.join(output_data, inp3)
            out=run([sys.executable,'script/boxplot_all_tumor_giusto.py',gene,feature,dir],shell=False, stdout=PIPE)
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
                    'dir':inp3,
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


####### DESEQ2 analisi###############

def choseimage(feature,tumor):
    pathfiles='rolls/static/media/deseq2/'+feature+'/'+tumor+'/'
    files=os.listdir(pathfiles)
    filelist=[]
    for file in files:
        if tumor in file:
            path='/media/deseq2/'+feature+'/'+tumor+'/'+file
            if 'jpg' in file:
                if 'EnhancedVolcano' in file:
                    enhancedimage=path
                else:
                    filelist.append(path)
        
    if len(filelist)>0:
        return(enhancedimage, filelist)
    else:
        return()



def deseq2copy1(request):
    if request.method == 'POST':
        form = Deseq2form(request.POST)
        if form.is_valid():
            
            tumor=form.cleaned_data['tumor']
            feature=form.cleaned_data['feature']
            
            dir='rolls/static/media/deseq2/'+feature+'/'+tumor
            if os.path.isdir(dir): 
                images=choseimage(feature,tumor)
                
                form=Deseq2form()              
                           
                return render(request, 'rolls/deseq2.html', {'form':form, 
                'feature': feature,
                'tumor':tumor,
                'enhancedimage': images[0],
                'images1': images[1][0],
                'images2': images[1][1],
                'images3': images[1][2],
                'go':'Valid',
                'parametri': parametri[feature],
                'dir':"http://160.80.35.91:7000/static/media/deseq2/"+feature+"/"+tumor,
                })
            else:
                form=Deseq2form()
                return render(request, 'rolls/deseq2.html', {'form':form,
                'feature': feature,
                'tumor':tumor, 
                'go':'error'})
        
    
    form=Deseq2form()
    return render(request, 'rolls/deseq2.html', {'form':form})


###################################################
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

            dir='rolls/static/media/deseq2/'+feature+'/'+tumor
            if os.path.isdir(dir): 
                images=choseimage(feature,tumor)
                

                form=Deseq2form()                         
                return render(request, 'rolls/deseq2.html', {'form':form, 
                    'feature': feature,
                    'tumor':tumor,
                    'enhancedimage': images[0],
                    'images1': images[1][0],
                    'images2': images[1][1],
                    'images3': images[1][2],
                    'go':'Valid',
                    'parametri': parametri[feature],
                    'dir':"http://160.80.35.91:7000/static/media/deseq2/"+feature+"/"+tumor,
                    })

            else:
                form=Deseq2form()
                return render(request, 'rolls/deseq2.html', {'form':form,
                'feature': feature,
                'tumor':tumor, 
                'go':'error'})



    form = Deseq2form()       
    return render(request, 'rolls/deseq2.html', {'form':form})





################################
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

