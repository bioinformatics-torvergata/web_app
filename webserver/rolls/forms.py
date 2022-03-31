from django import forms


class Analisiform(forms.Form):
    gene=forms.CharField()
    
    feature=forms.ChoiceField(
        choices=[
            ('none','Choice..'),
            ('gender','Gender'),
            ('age','Age'),
            ('radiation_therapy','Radiation therapy'),
            ('patient_status','Patient status'),
            ('diabetes', 'Diabetes'),
            ('tobacco_smoking_history','Tobacco smoking history'),
            ('menopause_status','Menopause status'),
            ('alcohol_history_documented','Alcohol history documented'),
            ('pathologic_stage','Pathologic stage'),
            ])
    tumor=forms.ChoiceField(
       
        choices=[
            
           ('none','Choice..'),
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
           ('UVM','UVM')
             ])
    

class Deseq2form(forms.Form):    
    feature=forms.ChoiceField(
        choices=[
            ('none','Choice..'),
            ('gender','Gender'),
            ('age_at_initial_pathologic_diagnosis','Age'),
            ('radiation_therapy','Radiation therapy'),
            ('patient_status','Patient status'),
            ('diabetes', 'Diabetes'),
            ('tobacco_smoking_history','Tobacco smoking history'),
            ('menopause_status','Menopause status'),
            ('alcohol_history_documented','Alcohol history documented'),
            ('pathologic_stage','Pathologic stage'),
            ])

    tumor=forms.ChoiceField(
        choices=[
           ('none','Choice..'),
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
           ('UVM','UVM')
             ])

    
class Analisiform1(forms.Form):
    gene=forms.CharField()
    tumor=forms.ChoiceField(
       
        choices=[
            
           ('none','Choice..'),
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
           ('UVM','UVM')
             ])
    
