from django.db import models

# Create your models here.
     
class Analisi(models.Model):
    gene= models.CharField(max_length=50)
    tumori=models.CharField(max_length=50)
    features=models.CharField(max_length=50)

