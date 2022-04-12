from django.urls import path
from . import views



urlpatterns=[
    path('', views.rolls, name='home'),
    path('documentation', views.documentation, name="documentation"),
    path('analisiprova', views.analisiprova, name= 'analisiprova'), 
    path('deseq2', views.deseq2, name='deseq2'),
    path('table', views.table, name='table'),
    
    path('overall_survival', views.overall_survival, name= 'overall_survival'),
    path('OS_pathway', views.os_pathway, name= 'OS_pathway'),
    path('OS_interaction', views.os_interaction, name='OS_interaction'),
    path('differential_expression', views.differential_expression, name= 'differential_expression'),
    
    
    path('downloadfile', views.downloadfile,name='downloadfile'),
    path('downloadfileszip', views.downloadfileszip,name='downloadfileszip')
    #url(r'^downloadfileszip/$', views.downloadfileszip,name='downloadfileszip'),

]
