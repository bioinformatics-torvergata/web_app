from django.shortcuts import render
from django.http import HttpResponse
from subprocess import run,PIPE
import sys
from matplotlib import image
from rolls.forms import Analisiform, Deseq2form, Analisiform1, Analisiformprova
import os
from django.http import StreamingHttpResponse
from wsgiref.util import FileWrapper
import mimetypes
from shutil import make_archive
import time

# Create your views here.


def rolls(request):
    return render(request, 'rolls/home.html')

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


####### analisi espressione differenziale #############
def differential_expression(request):
    if request.method == 'POST':
        form = Analisiformprova(request.POST)
        if form.is_valid():
            gene=form.cleaned_data['gene']
            tumor=form.cleaned_data['tumor']
            
            print(gene, tumor)
            inp1=gene
            inp2=tumor
            inp3=(time.strftime("%H%M%S%m"))
            out=run([sys.executable,'script/boxplot_all_tumor_giusto.py',inp1,inp2,inp3],shell=False, stdout=PIPE)
            print(out)
            
            dir='rolls/static/media/saveanalisi/'+inp3+'/'
            files=os.listdir(dir)
            for file in files:
                if file[-3:]=='png':
                    image='/media/saveanalisi/'+inp3+'/'+file
            form=Analisiformprova()
            return render(request, 'rolls/differential_expression.html', {
                'form':form, 
                'formresult': out.stdout.decode('ascii'),
                'image': image,
                'go':True,})

    form=Analisiformprova()
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
            
            dir='rolls/static/media/saveanalisi/overall_survival/'+inp3+'/'
            files=os.listdir(dir)
            for file in files:
                if file[-3:]=='png':
                    image='/media/saveanalisi/overall_survival/'+inp3+'/'+file
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
    pathfiles='/home/chiara/webserver/rolls/static/media/deseq2/'+feature+'/result/'+tumor+'/'
    files=os.listdir(pathfiles)
    filelist=[]
    for file in files:
        if tumor in file:
            path='/media/deseq2/'+feature+'/result/'+tumor+'/'+file
            if 'jpg' in file:
                if 'EnhancedVolcano' in file:
                    enhancedimage=path
                else:
                    filelist.append(path)
        
    if len(filelist)>0:
        return(enhancedimage, filelist)
    else:
        return()


def deseq2(request):
    if request.method == 'POST':
        form = Deseq2form(request.POST)
        if form.is_valid():
            tumor=form.cleaned_data['tumor']
            feature=form.cleaned_data['feature']

            images=choseimage(feature,tumor)
            
            form=Deseq2form()
            print(images)
            
            return render(request, 'rolls/deseq2.html', {'form':form, 
            'feature': feature,
            'tumor':tumor,
            'enhancedimage': images[0],
            'images1': images[1][0],
            'images2': images[1][1],
            'images3': images[1][2],
            'go':True,
            })

    form=Deseq2form()
    return render(request, 'rolls/deseq2.html', {'form':form})



#####download file zip ################
def downloadfileszip(request, file_name=''):
    """
    A django view to zip files in directory and send it as downloadable response to the browser.
    Args:
      @request: Django request object
      @file_name: Name of the directory to be zipped
    Returns:
      A downloadable Http response
    """
    files_path = "rolls/static/media/deseq2/gender/result/KIRP"
    path_to_zip = make_archive(files_path, "zip", files_path)
    response = HttpResponse(FileWrapper(open(path_to_zip,'rb')), content_type='application/zip')
    response['Content-Disposition'] = 'attachment; filename="{filename}.zip"'.format(
        filename = file_name.replace(" ", "_")
    )
    return response



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

