from django.urls import path
from . import views



urlpatterns=[
    path('', views.rolls, name='home'),
    path('documentation', views.documentation, name="documentation"),
    path('dataset', views.dataset, name="dataset"),
    path('contact', views.contact, name="contact"),
    path('correlation_analysis', views.correlation_analysis, name= 'correlation_analysis'), 
    path('deseq2', views.deseq2, name='deseq2'),
       
    path('overall_survival', views.overall_survival, name= 'overall_survival'),
    path('OS_pathway', views.os_pathway, name= 'OS_pathway'),
    path ('pathwayPROVA', views.pathwayPROVA, name='pathwayPROVA'),

    path('OS_interaction', views.os_interaction, name='OS_interaction'),
    
    
    
    path('differential_expression', views.differential_expression, name= 'differential_expression'),
    path('diff_exp_single_tumor', views.diff_exp_single_tumor, name= 'diff_exp_single_tumor'),
    path('rna_interactor',views.correlation_analysis, name= 'rna_interactor'),
    path('gene-suggestions/', views.gene_suggestions, name='gene_suggestions'),
    

]
