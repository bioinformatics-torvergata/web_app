from django.db import models

# Create your models here.

TUMOR=[
    (None,'Choice..'),
    ('ACC','ACC'),
    ('BLCA','BLCA'),
    ('BRCA','BRCA'),
    ('CESC','CESC'),
    ('CHOL','CHOL'),
    ('COAD','COAD'),
    ('DLBC','DLBC'),
    ('ESCA','ESCA'),
    ('GBM','GBM'),
    ('HNSC','HNSC'),
    ('KICH','KICH'),
    ('KIRC','KIRC'),
    ('KIRP','KIRP'),
    ('LGG','LGG'),
    ('LIHC','LIHC'),
    ('LUAD','LUAD'),
    ('LUSC','LUSC'),
    ('MESO','MESO'),
    ('OV','OV'),
    ('PAAD','PAAD'),
    ('PCPG','PCPG'),
    ('PRAD','PRAD'),
    ('READ','READ'),
    ('SARC','SARC'),
    ('SKCM','SKCM'),
    ('STAD','STAD'),
    ('TGCT','TGCT'),
    ('THCA','THCA'),
    ('THYM','THYM'),
    ('UCEC','UCEC'),
    ('UCS','UCS'),
    ('UVM','UVM')]

#rimosso MESO di cui non abbiamo i dati mutazionali
TUMOR_MUTATION=[ 
    (None,'Choice..'),
    ('ACC','ACC'),
    ('BLCA','BLCA'),
    ('BRCA','BRCA'),
    ('CESC','CESC'),
    ('CHOL','CHOL'),
    ('COAD','COAD'),
    ('DLBC','DLBC'),
    ('ESCA','ESCA'),
    ('GBM','GBM'),
    ('HNSC','HNSC'),
    ('KICH','KICH'),
    ('KIRC','KIRC'),
    ('KIRP','KIRP'),
    ('LGG','LGG'),
    ('LIHC','LIHC'),
    ('LUAD','LUAD'),
    ('LUSC','LUSC'),
    ('OV','OV'),
    ('PAAD','PAAD'),
    ('PCPG','PCPG'),
    ('PRAD','PRAD'),
    ('READ','READ'),
    ('SARC','SARC'),
    ('SKCM','SKCM'),
    ('STAD','STAD'),
    ('TGCT','TGCT'),
    ('THCA','THCA'),
    ('THYM','THYM'),
    ('UCEC','UCEC'),
    ('UCS','UCS'),
    ('UVM','UVM')]

FEATURES=[
    (None,'Choice..'),
    ('gender','Gender'),
    ('age_at_initial_pathologic_diagnosis','Age'),
    ('radiation_therapy','Radiation therapy'),
    ('patient_status','Patient status'),
    ('diabetes', 'Diabetes'),
    ('tobacco_smoking_history','Tobacco smoking history'),
    ('menopause_status','Menopause status'),
    ('alcohol_history_documented','Alcohol history documented'),
    ('pathologic_stage','Pathologic stage'),
    ]

FEATURE_R=[
    (None,'Please select a tumor type first'),
    ('gender','Gender'),
    ('alcohol_history_documented','Alcohol history documented'),
    # ('history_of_diabetes', 'Diabetes'),
    ('radiation_therapy','Radiation therapy'),
    ('person_neoplasm_cancer_status','Neoplasm cancer status'),
    ]

CHOICE_FEATURE =[
    ('Menopause status',(
        ('BRCA','BRCA'),
        ('OV','OV'),
        )
    )
]



class Analisi(models.Model):
    gene= models.CharField(
        max_length=20,)
    
    miRNA= models.CharField(
        max_length=20,)

    tumor=models.CharField(
        max_length=10,
        choices= TUMOR,
        default='Choice..',)

    feature=models.CharField(
        max_length=50,
        choices= FEATURES,
        default= 'Choice..',)
        

class Analisi_mutation(models.Model):
    gene= models.CharField(
        max_length=20,)
    
    tumor=models.CharField(
        max_length=10,
        choices= TUMOR_MUTATION,
        default='Choice..',)

    feature=models.CharField(
        max_length=50,
        choices= FEATURE_R,
        default= 'Choice..',) 

    NUMBER_CHOICES = [
        (10, '10'),
        (15, '15'),
        (20, '20'),
        (25, '25'),
    ]
    number = models.IntegerField(choices=NUMBER_CHOICES, default=10)

    feature_selected=models.CharField(
        max_length=50,
        default= 'Choice..',)

#class Downloadfileszip(models.Model)

#class link(models.Model):
#    url=models.TextField()


class Gene(models.Model):
    gene= models.CharField(max_length=100)


class Pathway(models.Model):
    pathway= models.CharField(max_length=100)

class Protein(models.Model):
    protein= models.CharField(max_length=100)

class Gene_symbol(models.Model):
    gene_symbol= models.CharField(max_length=100)