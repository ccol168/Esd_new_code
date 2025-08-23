import os
import re
import argparse
from collections import defaultdict

def slice_run (filename) :

    slicelength = 100
    store_lines = []
    with open(filename) as file :
        for line in file :
            store_lines.append(line)
    j = 0
    for i in range(0, len(store_lines), slicelength):
        chunk = store_lines[i:i + slicelength]

        outfilename = filename[:-5] + "_" + str(j) + ".list"
    
        with open(outfilename, "w") as fileout:
            for item in chunk:
                fileout.write(str(item))  # Write each element on a new line

        j = j + 1

def processa_files(cartella_completa, cartella_output, ds_tag, cartella_finale,end_file):
    """Processa i file con il tag specificato e scrive i risultati nelle liste."""
    
    pattern_run = re.compile(r"RUN\.(\d+)\.")
    run_files = defaultdict(list)

    for root, _, files in os.walk(cartella_completa):
        files.sort()  # Ordina i file per nome
        for file in files:
            if file.endswith(end_file) and ds_tag in file:
                match = pattern_run.search(file)
                if match:
                    numero_run = match.group(1)
                    percorso_completo = os.path.join(root, file)
                    run_files[numero_run].append(percorso_completo)

    if not run_files:
        print(f"Nessun file trovato per {ds_tag} in {cartella_completa}.")
        return

    os.makedirs(cartella_output, exist_ok=True)  # Crea la cartella se non esiste

    for numero_run, files in run_files.items():
        nome_file_output = os.path.join(cartella_output, f"rawfile_RUN{numero_run}_2025{cartella_finale}.list")
        try:
            with open(nome_file_output, 'w') as output_file:
                output_file.write('\n'.join(files) + '\n')
            print(f"Creato file: {nome_file_output} con {len(files)} file")
            slice_run(nome_file_output)
            os.remove(nome_file_output)
        except IOError as e:
            print(f"Errore nella scrittura del file {nome_file_output}: {e}")

def crea_liste_per_run(cartella_finale,junosw_version):
    """Crea file di lista per RUN, processando sia 'ds-2' che 'ds-3'."""
    
    percorso_base = f"/storage/gpfs_data/juno/junofs/production/storm/dirac/juno/esd/{junosw_version}/2025/"
    cartella_completa = os.path.join(percorso_base, cartella_finale)

    percorso_rtraw = "/storage/gpfs_data/juno/junofs/production/storm/dirac/juno/rtraw/2025/"
    cartella_rtraw = os.path.join (percorso_rtraw,cartella_finale)

    if not os.path.exists(cartella_completa):
        print(f"Errore: La cartella {cartella_completa} non esiste.")
        return

    cartella_list_ds2 = os.path.join(os.getcwd(), "list")
    cartella_list_ds3 = os.path.join(os.getcwd(), "list_mm")
    cartella_list_rtraw = os.path.join(os.getcwd(), "list_rtraw")

    # Processa entrambi i tipi di file includendo il nome della sottocartella
    processa_files(cartella_completa, cartella_list_ds2, "ds-2", cartella_finale,".esd")
    processa_files(cartella_completa, cartella_list_ds3, "ds-3", cartella_finale,".esd")
    processa_files(cartella_rtraw, cartella_list_rtraw, "ds-2", cartella_finale,".rtraw")

def main():
    parser = argparse.ArgumentParser(description="Crea file .list per RUN dai file .rtraw.")
    parser.add_argument("cartella_finale", help="Nome della sottocartella finale (es. '1208').")
    parser.add_argument("-junosw_version",help="junosw version to find the esd files",default="J25.4.1.b")

    args = parser.parse_args()
    crea_liste_per_run(args.cartella_finale,args.junosw_version)

if __name__ == "__main__":
    main()

