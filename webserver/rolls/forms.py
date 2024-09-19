from django import forms
from .models import Analisi, Analisi_mutation

class Gene(forms.ModelForm):
    gene=forms.CharField(widget=forms.TextInput(attrs={"placeHolder":"gene symbol/ENSG"}))
    class Meta:
        model=Analisi
        fields=('gene','tumor')

class Analisiformcompleto(forms.ModelForm):
    #gene=forms.CharField(widget=forms.TextInput(attrs={"placeHolder":"gene(ENSG)/miRNA/protein"}))
    gene = forms.CharField(
        widget=forms.TextInput(
            attrs={
                "placeholder": "Gene(ENSG)/miRNA/protein",
                "id": "gene-input",  # Aggiungi l'id qui per collegarlo all'autocomplete
                "autocomplete": "on",  # Disabilita l'autocomplete del browser
            }
        )
    )
    class Meta:
        model= Analisi
        fields=('gene', 'tumor','feature')


class Analisiform(forms.ModelForm):
    gene = forms.CharField(
        widget=forms.TextInput(
            attrs={
                "placeholder": "Gene(ENSG)/miRNA/protein",
                "id": "gene-input",  # Aggiungi l'id qui per collegarlo all'autocomplete
                "autocomplete": "on",  # Disabilita l'autocomplete del browser
            }
        )
    )
    class Meta:
        model= Analisi
        fields=('gene','feature')




class Analisiform1(forms.ModelForm):
    gene = forms.CharField(
        widget=forms.TextInput(
            attrs={
                "placeholder": "Gene(ENSG)/miRNA/protein",
                "id": "gene-input",  # Aggiungi l'id qui per collegarlo all'autocomplete
                "autocomplete": "on",  # Disabilita l'autocomplete del browser
            }
        )
    )

    class Meta:
        model = Analisi
        fields = ('gene', 'tumor')

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
    
    pathway=forms.CharField(
        widget=forms.TextInput(
            attrs={
                "placeHolder":"es. SA_PROGRAMMED_CELL_DEATH ",
                "id": "pathway-input",  # Aggiungi l'id qui per collegarlo all'autocomplete
                "autocomplete": "on",  # Disabilita l'autocomplete del browser
                }))
    class Meta:
        model= Analisi
        fields=('tumor',)   


class FormTumorMutation(forms.ModelForm): 

    class Meta:
        model=Analisi_mutation
        fields=('tumor',)

class Onlygeneform(forms.ModelForm):
    gene = forms.CharField(
        widget=forms.TextInput(
            attrs={
                "placeholder": "Gene (Gene Symbol)",
                "id": "gene-input",  # Aggiungi l'id qui per collegarlo all'autocomplete
                "autocomplete": "on",  # Disabilita l'autocomplete del browser
            }
        )
    )
    class Meta:
        model= Analisi
        fields=('gene',)

