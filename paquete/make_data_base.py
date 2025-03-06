import os
import requests
from paquete import process_pdb
from paquete import calculate_angles
from paquete import get_atom_coords
from paquete import torsion_angle
import pandas as pd


def create_pdb_folder(save_path="pdb_files/"):
    """Crea la carpeta donde se guardarán los archivos PDB si no existe."""
    if not os.path.exists(save_path):
        os.makedirs(save_path)
        print(f"Carpeta creada: {save_path}")
    else:
        print(f"La carpeta {save_path} ya existe.")


def download_pdb(pdb_id, save_path="pdb_files/"):
    """
    Descarga un archivo PDB dado su identificador y lo guarda en la carpeta especificada.

    Args:
    - pdb_id (str): Identificador del archivo PDB (Ej: '1CRN').
    - save_path (str): Carpeta donde se guardarán los archivos descargados.
    """
    create_pdb_folder(save_path)  # Asegurar que la carpeta existe

    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    response = requests.get(url)

    if response.status_code == 200:
        with open(f"{save_path}{pdb_id}.pdb", "w") as file:
            file.write(response.text)
        print(f"PDB {pdb_id} descargado exitosamente.")
    else:
        print(f"Error al descargar {pdb_id}. Verifica si el ID es correcto.")


def download_pdbs_from_txt(txt_file, save_path="pdb_files/"):
    """
    Lee un archivo .txt con una lista de IDs de PDB y los descarga automáticamente.

    Args:
    - txt_file (str): Ruta al archivo .txt con los nombres de los PDB.
    - save_path (str): Carpeta donde se guardarán los archivos descargados.
    """
    create_pdb_folder(save_path)  # Crear la carpeta antes de descargar

    if not os.path.exists(txt_file):
        print("Archivo de PDBs no encontrado.")
        return

    with open(txt_file, "r") as file:
        pdb_ids = [line.strip() for line in file.readlines() if line.strip()]

    print(f"Descargando {len(pdb_ids)} archivos PDB...")

    for pdb_id in pdb_ids:
        download_pdb(pdb_id, save_path)

    print("Descarga completada.")


# download_pdbs_from_txt("listapdb.txt")


def process_multiple_pdbs(pdb_folder):
    """
    Procesa múltiples archivos PDB para extraer los ángulos phi y psi.

    Args:
    - pdb_folder (str): Carpeta donde se encuentran los archivos PDB.

    Returns:
    - pd.DataFrame: DataFrame con los ángulos phi y psi de todos los PDB procesados.
    """
    all_angles = []

    for pdb_file in os.listdir(pdb_folder):
        if pdb_file.endswith(".pdb"):
            file_path = os.path.join(pdb_folder, pdb_file)
            df = process_pdb(file_path)  # Extraer información del PDB
            angles_df = calculate_angles(df, get_atom_coords, torsion_angle)  # Calcular ángulos
            angles_df["pdb_id"] = pdb_file  # Añadir el nombre del PDB
            all_angles.append(angles_df)

    return pd.concat(all_angles, ignore_index=True)


# Aplicar la función a la carpeta "pdb_files/"
angles_data = process_multiple_pdbs("pdb_files/")
angles_data.to_csv("phi_psi_angles_dataset.csv", index=False)

# print("Archivo con ángulos phi y psi guardado: phi_psi_angles_dataset.csv")