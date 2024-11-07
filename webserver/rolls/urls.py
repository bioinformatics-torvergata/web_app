from django.urls import path
from . import views



urlpatterns=[
    path('', views.rolls, name='home'),
    
    path('documentation', views.documentation, name="documentation"),
    path('dataset', views.dataset, name="dataset"),
    path('contact', views.contact, name="contact"),
    path('correlation_analysis', views.correlation_analysis, name= 'correlation_analysis'), 
   
    path('deseq2', views.deseq2, name='deseq2'),
    path('differential_expression', views.differential_expression, name= 'differential_expression'),
    path('diff_exp_single_tumor', views.diff_exp_single_tumor, name= 'diff_exp_single_tumor'),
    path('get-features/', views.get_features, name='get_features'),  # L'URL per la richiesta AJAX
    path('gene-suggestions/', views.gene_suggestions, name='gene_suggestions'),

    path('differential_expression_protein', views.differential_expression_protein, name= 'differential_expression_protein'),
    path('diff_exp_single_tumor_protein', views.diff_exp_single_tumor_protein, name= 'diff_exp_single_tumor_protein'),
    path('protein_suggestions/', views.protein_suggestions, name='protein_suggestions'),
       
    path('overall_survival', views.overall_survival, name= 'overall_survival'),
    path('OS_pathway', views.os_pathway, name= 'OS_pathway'),
    path('pathwayPROVA', views.pathwayPROVA, name='pathwayPROVA'),
    path('pathway_suggestions/', views.pathway_suggestions, name='pathway_suggestions'),

    path('rna_interactor',views.correlation_analysis, name= 'rna_interactor'),
    path('OS_interaction', views.os_interaction, name='OS_interaction'),
    
    path('tumor_mutation_analysis/', views.tumor_mutation_analysis, name='tumor_mutation_analysis'),
    path('tumor_oncoplot/', views.tumor_oncoplot, name='tumor_oncoplot'),
    path('somatic_interaction_analysis/', views.somatic_interaction_analysis, name='somatic_interaction_analysis'),
    path('gene_mutation_analysis/', views.gene_mutation_analysis, name='gene_mutation_analysis'),
    path('de_mut/', views.de_mut, name='de_mut'),
    path('de_mut_clinical_feature/', views.de_mut_clinical_feature, name='de_mut_clinical_feature'),
    path('gene_symbol_suggestions/', views.gene_symbol_suggestions, name='gene_symbol_suggestions'),
    path('survival_with_gene_mutation_status/', views.survival_with_gene_mutation_status, name='survival_with_gene_mutation_status'),
    path('get-features_R/', views.get_features_R, name='get_features_R'),  # L'URL per la richiesta AJAX


    path('deconvolution/', views.deconvolution,name='deconvolution'),
    path('corr_cell_pathway/', views.corr_cell_pathway,name='corr_cell_pathway'),

    #path('download_zip/', views.download_zip, name='download_zip'),
    path('download_zip/<str:subdirectory>/<str:zip_name>/', views.download_zip, name='download_zip'),

]
