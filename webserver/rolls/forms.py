from django import forms
from .models import Analisi


class Analisiformprova(forms.ModelForm):
    gene=forms.CharField(widget=forms.TextInput(attrs={"placeHolder":"gene/miRNA/protein"}))
    class Meta:
        model= Analisi
        fields=('gene', 'tumor','feature')


class Analisiform(forms.ModelForm):
    gene=forms.CharField(widget=forms.TextInput(attrs={"placeHolder":"gene/miRNA/protein"}))
    class Meta:
        model= Analisi
        fields=('gene','feature')



class Analisiform1(forms.ModelForm):
    gene=forms.CharField(widget=forms.TextInput(attrs={"placeHolder":"gene/miRNA/protein"}))

    class Meta:
        model= Analisi
        fields=('gene', 'tumor')   


class Deseq2form(forms.ModelForm): 

    class Meta:
        model=Analisi
        fields=('tumor','feature')

   
class AnalisiOSinteraction(forms.ModelForm):
    gene=forms.CharField(widget=forms.TextInput(attrs={"placeHolder":"ENSG/gene symbol"}))
    class Meta:
        model=Analisi
        fields=('gene','miRNA', 'tumor')


   
