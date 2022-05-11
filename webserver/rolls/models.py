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
    ('LUSH','LUSH'),
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
        
        
        
#class Downloadfileszip(models.Model)

#class link(models.Model):
#    url=models.TextField()
