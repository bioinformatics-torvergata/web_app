from django.shortcuts import render
from django.http import HttpResponse
from subprocess import run,PIPE
import sys
from matplotlib import image
from rolls.forms import Analisiform, Deseq2form, Analisiform1, Analisiformprova, AnalisiOSinteraction
import os
from django.http import StreamingHttpResponse
from wsgiref.util import FileWrapper
import mimetypes
from shutil import make_archive
import time
import os.path



def rolls(request):
    return render(request, 'rolls/home.html')

def documentation(request):
    return render(request, 'rolls/documentation.html')

def table(request):
    return render(request, 'rolls/table.html')


def os_pathway(request):
    return render(request, 'rolls/OS_pathway.html')

#############Overall survival con dati interazione miRNA-mRNA ############
def os_interaction(request):
    if request.method == 'POST':
        form = AnalisiOSinteraction(request.POST)
        if form.is_valid():
            gene=form.cleaned_data['gene']
            miRNA=form.cleaned_data['miRNA']
            tumor= form.cleaned_data['tumor']

            inp3=(time.strftime("%H%M%S%m"))
            out=run([sys.executable,'script/overall_survival_interaction.py',gene,miRNA,tumor,inp3],shell=False, stdout=PIPE)
            print(out)
            dir='rolls/static/media/saveanalisi/'+inp3+'/'
            if os.path.isdir(dir): 
                #dir='rolls/static/media/saveanalisi/'+inp3+'/'
                files=os.listdir(dir)
                images=[]
                for file in files:
                    if file[-3:]=='jpg':
                        images.append('/media/saveanalisi/'+inp3+'/'+file)

                form=AnalisiOSinteraction()
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
                form=AnalisiOSinteraction()
                return render(request, 'rolls/OS_interaction.html', {'form':form,
                'miRNA':miRNA,
                'gene':gene,
                'tumor':tumor, 
                'go':'error'})
                

    form=AnalisiOSinteraction()
    return render(request, 'rolls/OS_interaction.html', {'form':form})
    


def analisiprova(request):
    if request.method == 'POST':
        form = Analisiform(request.POST)
        if form.is_valid():
            gene=form.cleaned_data['gene']
            tumor=form.cleaned_data['tumor']
            feature=form.cleaned_data['feature']
            print(gene, tumor, feature)
            inp1=gene
            inp2=tumor
            inp3=(time.strftime("%H%M%S%m"))
            out=run([sys.executable,'script/#',inp1,inp2,inp3],shell=False, stdout=PIPE)
            print(out)
            
            dir='rolls/static/media/saveanalisi/'+inp3+'/'
            files=os.listdir(dir)
            for file in files:
                if file[-3:]=='png':
                    image='/media/saveanalisi/'+inp3+'/'+file
            form=Analisiform()
            return render(request, 'rolls/form.html', {
                'form':form, 
                'formresult': out.stdout.decode('ascii'),
                'image': image,
                'go':True,})

    form=Analisiform()
    return render(request, 'rolls/form.html', {'form':form})


####### differential expression analisys with all tumor for feature #############
def differential_expression(request):
    if request.method == 'POST':
        form = Analisiform(request.POST)
        if form.is_valid():
            gene=form.cleaned_data['gene']
            feature=form.cleaned_data['feature']
            
            inp3=(time.strftime("%H%M%S%m"))
            out=run([sys.executable,'script/boxplot_all_tumor_giusto.py',gene,feature,inp3],shell=False, stdout=PIPE)
            print(out)
            dir='rolls/static/media/saveanalisi/'+inp3+'/'
            files=os.listdir(dir)
            for file in files:
                if file[-3:]=='jpg':
                    image='/media/saveanalisi/'+inp3+'/'+file
            form=Analisiform()
            return render(request, 'rolls/differential_expression.html', {
                'form':form, 
                'formresult': out.stdout.decode('ascii'),
                'image': image,
                'go':True,
                'gene':gene,
                'feature':feature,})

    form=Analisiform()
    return render(request, 'rolls/differential_expression.html', {'form':form})



######### OVERALL SURVIVAL ################
def overall_survival(request):
    if request.method == 'POST':
        form = Analisiform1(request.POST)
        if form.is_valid():
            gene=form.cleaned_data['gene']
            tumor=form.cleaned_data['tumor']
            print(gene, tumor)
            
            inp3=(time.strftime("%H%M%S%m"))
            out=run([sys.executable,'script/overall_survival.py',gene,tumor,inp3],shell=False, stdout=PIPE)
            print(out)
            
            dir='rolls/static/media/saveanalisi/'+inp3+'/'
            files=os.listdir(dir)
            for file in files:
                if file[-3:]=='png':
                    image='/media/saveanalisi/'+inp3+'/'+file
            form=Analisiform1()
            return render(request, 'rolls/overall_survival.html', {
                'form':form, 
                'formresult': out.stdout.decode('ascii'),
                'image': image,
                'go':True,
                'gene':gene ,
                'tumor':tumor,})

    form=Analisiform1()
    return render(request, 'rolls/overall_survival.html', {'form':form})


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

parametri={
    'patient_status':'Tumor vs Ctrl',
    'gender':'Female vs Male',
    'pathologic_stage': 'StageIII_IV vs StageI_II',
    'radiation_therapy': 'YES vs NO',
    'diabetes': 'YES vs NO',
    'tobacco_smoking_history': 'Smoker vs Non_Smoker',
    'menopause_status':'Post-menopause vs Pre-menopause',
    'alcohol_history_documented': 'YES vs NO',
    'age_at_initial_pathologic_diagnosis': 'Above the median vs Below the median' ,
}

def deseq2(request):
    if request.method == 'POST':
        form = Deseq2form(request.POST)
        if form.is_valid():
            
            tumor=form.cleaned_data['tumor']
            feature=form.cleaned_data['feature']
            
            cartella='rolls/static/media/deseq2/'+feature+'/'+tumor
            if os.path.isdir(cartella): 
                images=choseimage(feature,tumor)
                
                form=Deseq2form()
                print(images)
                
                downloadfileszip(tumor)
                
                
                return render(request, 'rolls/deseq2.html', {'form':form, 
                'feature': feature,
                'tumor':tumor,
                'enhancedimage': images[0],
                'images1': images[1][0],
                'images2': images[1][1],
                'images3': images[1][2],
                'go':'Valid',
                'parametri': parametri[feature],
                })
            else:
                form=Deseq2form()
                return render(request, 'rolls/deseq2.html', {'form':form,
                'feature': feature,
                'tumor':tumor, 
                'go':'error'})
        
    
    form=Deseq2form()
    return render(request, 'rolls/deseq2.html', {'form':form})



#cartella='rolls/static/media/deseq2/gender/'+tumor
#####download file zip ################
def downloadfileszip(request,tumor):
    file_name=tumor
    files_path = 'rolls/static/media/deseq2/gender/'+tumor
    path_to_zip = make_archive(files_path, "zip", files_path)
    response = HttpResponse(FileWrapper(open(path_to_zip,'rb')), content_type='application/zip')
    response['Content-Disposition'] = 'attachment; filename="{filename}.zip"'.format(
        filename = file_name.replace(" ", "_")
    )
    return response




'''
def downloadfileszip(request, file_name=''):
    """
    A django view to zip files in directory and send it as downloadable response to the browser.
    Args:
      @request: Django request object
      @file_name: Name of the directory to be zipped
    Returns:
      A downloadable Http response
    """
    #file_name=tumor
    files_path = 'rolls/static/media/deseq2/gender/'+tumor
    path_to_zip = make_archive(files_path, "zip", files_path)
    response = HttpResponse(FileWrapper(open(path_to_zip,'rb')), content_type='application/zip')
    response['Content-Disposition'] = 'attachment; filename="{filename}.zip"'.format(
        filename = file_name.replace(" ", "_")
    )
    return response
'''


#download button
def downloadfile(request):
    base_dir=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    filename='EnhancedVolcanoGeni_KIRP.jpg'
    filepath=base_dir + '/files/'+ filename
    thefile=filepath
    filename=os.path.basename(thefile)
    chunk_size= 8192
    response= StreamingHttpResponse(FileWrapper(open(thefile,'rb'), chunk_size),content_type=mimetypes.guess_type(thefile)[0])
    response['Content-Lenght']=os.path.getsize(thefile)
    response['Content-Disposition']="Attachment;filename=%s" %filename
    return response

