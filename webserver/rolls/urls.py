from django.urls import path
from . import views



urlpatterns=[
    path('', views.rolls, name='home'),
    path('analisiprova', views.analisiprova, name= 'analisiprova'), 
    path('deseq2', views.deseq2, name='deseq2'),
    path('downloadfile', views.downloadfile,name='downloadfile'),
    path('downloadfileszip/dirname', views.downloadfileszip,name='downloadfileszip'),
    path('overall_survival', views.overall_survival, name= 'overall_survival')
]
