import os
import django

# Imposta il contesto di Django
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'webserver.settings')  # Sostituisci con il nome del tuo progetto
django.setup()

from rolls.models import Gene_symbol  # Sostituisci con il nome della tua app

# Funzione per popolare il modello Gene
def populate_genes():
    # Leggi la lista di geni dal file
    with open('gene_symbol_list.txt', 'r') as f:
        genes = f.readlines()

    # Rimuovi eventuali spazi bianchi e newline
    genes = [gene.strip() for gene in genes]

    # Inserisci ciascun gene nel database
    for gene_name in genes:
        gene_symbol, created = Gene_symbol.objects.get_or_create(gene_symbol=gene_name)
        if created:
            print(f'Creato gene: {gene_name}')
        else:
            print(f'Il gene {gene_name} esiste già.')

if __name__ == '__main__':
    populate_genes()
