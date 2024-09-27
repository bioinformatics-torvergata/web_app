import os
import django

# Imposta il contesto di Django
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'webserver.settings')  # Sostituisci con il nome del tuo progetto
django.setup()

from rolls.models import Protein  # Sostituisci con il nome della tua app

# Funzione per popolare il modello Gene
def populate_genes():
    # Leggi la lista di geni dal file
    with open('protein_list.txt', 'r') as f:
        genes = f.readlines()

    # Rimuovi eventuali spazi bianchi e newline
    proteins = [gene.strip() for gene in genes]

    # Inserisci ciascun gene nel database
    for p in proteins:
        protein, created = Protein.objects.get_or_create(protein=p)
        if created:
            print(f'Creato gene: {p}')
        else:
            print(f'Il gene {p} esiste già.')

if __name__ == '__main__':
    populate_genes()
