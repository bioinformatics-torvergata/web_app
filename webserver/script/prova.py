
import os
import configparser
from pathlib import Path

stringa='a,b,c'

#print(stringa)


# Carica il file di configurazione
config = configparser.ConfigParser()
script_dir = Path(__file__).parent

# Costruisci il percorso relativo al file di configurazione
config_file = script_dir.parent.parent / 'webserver' / 'webserver' / 'conf.ini'

config.read(config_file) 
#directory base
output_data = config['Paths']['output_data']
base_dir = config.get('Paths', 'base_dir', fallback='')

# Costruisci i percorsi completi
def get_full_path(relative_path):
    return os.path.join(base_dir, relative_path)

# Accesso ai valori e costruzione dei percorsi
if 'miRNA' in config:
    miRNA_name = get_full_path(config['miRNA']['name'])
    miRNA_dataframe = get_full_path(config['miRNA']['dataframe'])
    print(f"miRNA name path: {miRNA_name}")
    print(f"miRNA dataframe path: {miRNA_dataframe}")

if 'gene' in config:
    gene_name_ENSG = get_full_path(config['gene']['name_ENSG'])
    gene_dataframe = get_full_path(config['gene']['dataframe'])
    gene_dataframe_FPKM = get_full_path(config['gene']['dataframe_FPKM'])
    print(f"Gene ENSG path: {gene_name_ENSG}")
    print(f"Gene dataframe path: {gene_dataframe}")
    print(f"Gene dataframe FPKM path: {gene_dataframe_FPKM}")

if 'protein' in config:
    protein_name = get_full_path(config['protein']['name'])
    protein_dataframe = get_full_path(config['protein']['dataframe'])
    print(f"Protein name path: {protein_name}")
    print(f"Protein dataframe path: {protein_dataframe}")

if 'clinical' in config:
    clinical_data= get_full_path(config['clinical']['dati_clinici'])
    clinical_origin=get_full_path(config['clinical']['clinical_OS'])


if 'os' in config:
    os_pathway=get_full_path(config['OS']['os_pathway'])
    