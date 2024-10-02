import os
import django

# Set up the Django context
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'webserver.settings')  # Sostituisci con il nome del tuo progetto
django.setup()


from rolls.models import Gene_symbol  # Replace with your app name

# Function to delete genes from the model
def delete_genes_from_file(file_path):
    # Read the list of genes from the file
    with open(file_path, 'r') as f:
        genes_to_delete = f.readlines()

    # Remove any whitespace and newline characters
    genes_to_delete = [gene.strip() for gene in genes_to_delete]

    # Delete each gene from the database
    for gene_name in genes_to_delete:
        deleted, _ = Gene_symbol.objects.filter(gene_symbol=gene_name).delete()
        
        if deleted:
            print(f'Deleted gene: {gene_name}')
        else:
            print(f'The gene {gene_name} does not exist.')

if __name__ == '__main__':
    # Pass the file name as an argument
    file_name = 'gene_symbol_list.txt'  # Replace with your actual file name
    delete_genes_from_file(file_name)
    total_objects = Gene_symbol.objects.count()
    print(f"Numero totale di oggetti: {total_objects}")
