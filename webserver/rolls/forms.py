from django import forms
from .models import Analisi

class Gene(forms.ModelForm):
    gene=forms.CharField(widget=forms.TextInput(attrs={"placeHolder":"gene symbol/ENSG"}))
    class Meta:
        model=Analisi
        fields=('gene','tumor')

class Analisiformcompleto(forms.ModelForm):
    gene=forms.CharField(widget=forms.TextInput(attrs={"placeHolder":"gene(ENSG)/miRNA/protein"}))
    class Meta:
        model= Analisi
        fields=('gene', 'tumor','feature')


class Analisiform(forms.ModelForm):
    gene=forms.CharField(widget=forms.TextInput(attrs={"placeHolder":"gene(ENSG)/miRNA/protein"}))
    class Meta:
        model= Analisi
        fields=('gene','feature')



class Analisiform1(forms.ModelForm):
    gene=forms.CharField(widget=forms.TextInput(attrs={"placeHolder":"Gene(ENSG)/miRNA/protein"}))
    
    class Meta:
        model= Analisi
        fields=('gene','tumor')   


class Deseq2form(forms.ModelForm): 

    class Meta:
        model=Analisi
        fields=('tumor',)

   
class Analisi_interaction(forms.ModelForm):
    gene=forms.CharField(widget=forms.TextInput(attrs={"placeHolder":"ENSG/gene symbol"}))
    class Meta:
        model=Analisi
        fields=('gene','miRNA', 'tumor')


   
class Analisipath(forms.ModelForm):
    pathway=forms.CharField(widget=forms.TextInput(attrs={"placeHolder":"es. SA_PROGRAMMED_CELL_DEATH "}))
    class Meta:
        model= Analisi
        fields=('tumor',)   
