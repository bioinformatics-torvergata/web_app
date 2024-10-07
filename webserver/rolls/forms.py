from django import forms
from .models import Analisi, Analisi_mutation,TUMOR_MUTATION, FEATURES

class Gene(forms.ModelForm):
    gene=forms.CharField(widget=forms.TextInput(attrs={"placeHolder":"gene symbol/ENSG"}))
    class Meta:
        model=Analisi
        fields=('gene','tumor')

class Analisiformcompleto_old(forms.ModelForm):
    #gene=forms.CharField(widget=forms.TextInput(attrs={"placeHolder":"gene(ENSG)/miRNA/protein"}))
    gene = forms.CharField(
        widget=forms.TextInput(
            attrs={
                "placeholder": "ENSG/Gene symbol/miRNA",
                "id": "gene-input",  # Aggiungi l'id qui per collegarlo all'autocomplete
                "autocomplete": "on",  # Disabilita l'autocomplete del browser
            }
        )
    )
    class Meta:
        model= Analisi
        fields=('gene', 'tumor','feature')

###new
class Analisiformcompleto(forms.ModelForm):
    gene = forms.CharField(
        widget=forms.TextInput(
            attrs={
                "placeholder": "ENSG/Gene symbol/miRNA",
                "id": "gene-input",  # Aggiungi l'id qui per collegarlo all'autocomplete
                "autocomplete": "on",  # Disabilita l'autocomplete del browser
            }
        )
    )
    
    tumor = forms.ChoiceField(
        choices=TUMOR_MUTATION,
        widget=forms.Select(
            attrs={
                "id": "id_tumor"  # ID per il campo tumor, necessario per il JavaScript
            }
        )
    )
    
    feature = forms.ChoiceField(
        choices=FEATURES,
        widget=forms.Select(
            attrs={
                "id": "feature-select"  # ID per il campo feature
            }
        ),
        label="Select clinical feature"
    )
    
    class Meta:
        model = Analisi
        fields = ('gene', 'tumor', 'feature')

###############


class Analisiform(forms.ModelForm):
    gene = forms.CharField(
        widget=forms.TextInput(
            attrs={
                "placeholder": "ENSG/Gene symbol/miRNA",
                "id": "gene-input",  # Aggiungi l'id qui per collegarlo all'autocomplete
                "autocomplete": "on",  # Disabilita l'autocomplete del browser
            }
        )
    )
    class Meta:
        model= Analisi
        fields=('gene','feature')

#####FORM for PROTEOMICS ANALYSES #############

class Analisiform_protein(forms.ModelForm):
    protein = forms.CharField(
        widget=forms.TextInput(
            attrs={
                "placeholder": "insert protein",
                "id": "protein-input",  # Aggiungi l'id qui per collegarlo all'autocomplete
                "autocomplete": "on",  # Disabilita l'autocomplete del browser
            }
        )
    )
    class Meta:
        model= Analisi
        fields=('protein','feature')


class Analisiformcompleto_protein_old(forms.ModelForm):
    #gene=forms.CharField(widget=forms.TextInput(attrs={"placeHolder":"gene(ENSG)/miRNA/protein"}))
    protein = forms.CharField(
        widget=forms.TextInput(
            attrs={
                "placeholder": "insert protein",
                "id": "protein-input",  # Aggiungi l'id qui per collegarlo all'autocomplete
                "autocomplete": "on",  # Disabilita l'autocomplete del browser
            }
        )
    )
    class Meta:
        model= Analisi
        fields=('protein', 'tumor','feature')


class Analisiformcompleto_protein(forms.ModelForm):
    protein = forms.CharField(
        widget=forms.TextInput(
            attrs={
                "placeholder": "insert protein",
                "id": "protein-input",  # Aggiungi l'id qui per collegarlo all'autocomplete
                "autocomplete": "off",  # Disabilita l'autocomplete del browser
            }
        )
    )
    tumor = forms.ChoiceField(
        choices=TUMOR_MUTATION,
        widget=forms.Select(
            attrs={
                "id": "id_tumor"  # ID per il campo tumor, necessario per il JavaScript
            }
        )
    )
    
    feature = forms.ChoiceField(
        choices=FEATURES,
        widget=forms.Select(
            attrs={
                "id": "feature-select"  # ID per il campo feature
            }
        ),
        label="Select clinical feature"
    )
    class Meta:
        model= Analisi
        fields=('protein', 'tumor','feature')
#############################################################

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

class formSurvival(forms.ModelForm):
    gene = forms.CharField(
        widget=forms.TextInput(
            attrs={
                "placeholder": "Gene(ENSG)/miRNA/protein",
                "id": "gene-input",  # Aggiungi l'id qui per collegarlo all'autocomplete
                "autocomplete": "on",  # Disabilita l'autocomplete del browser
            }
        )
    )
    Methods = forms.ChoiceField(
        choices=[
            ('OS.time', 'OS.time'),
            ('DFI.time', 'DFI.time'),
        ],
         widget=forms.RadioSelect
    )


    class Meta:
        model = Analisi
        fields = ('gene', 'tumor','Methods')

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
    Methods = forms.ChoiceField(
        choices=[
            ('OS.time', 'OS.time'),
            ('DFI.time', 'DFI.time'),
        ],
         widget=forms.RadioSelect
    )
    class Meta:
        model= Analisi
        fields=('tumor','Methods')   


class FormTumorMutation(forms.ModelForm): 

    class Meta:
        model=Analisi_mutation
        fields=('tumor',)

class FormMutationChoice(forms.ModelForm):
    number = forms.ChoiceField(
        choices=[
            (10, '10'),
            (15, '15'),
            (20, '20'),
            (25, '25'),
        ],
        label="Select number of genes to be analysed"
    )

    class Meta:
        model = Analisi_mutation
        fields = ('tumor', 'number')

class tumorGeneform(forms.ModelForm):
    gene = forms.CharField(
        widget=forms.TextInput(
            attrs={
                "placeholder": "Gene (Gene Symbol)",
                "id": "gene_symbol-input",  # Aggiungi l'id qui per collegarlo all'autocomplete
                "autocomplete": "on",  # Disabilita l'autocomplete del browser
            }
        )
    )

    class Meta:
        model= Analisi_mutation
        fields=('gene','tumor')

