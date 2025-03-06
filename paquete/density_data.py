import pandas as pd
import numpy as np


# Cargar los datos de ángulos phi y psi
def load_phi_psi_data(csv_file):
    """Carga los datos de ángulos Ramachandran desde un archivo CSV."""
    return pd.read_csv(csv_file)


# Construcción de la tabla de densidad

def build_phi_psi_density_table(csv_file, bin_size=2):
    """Genera una tabla de densidad de ocupación para phi y psi sin clasificar por tipo de residuo."""
    angles_df = load_phi_psi_data(csv_file)

    # Discretizar los valores en bins de 2°
    angles_df["phi_bin"] = (angles_df["phi"] // bin_size) * bin_size
    angles_df["psi_bin"] = (angles_df["psi"] // bin_size) * bin_size

    # Contar ocurrencias
    density_table = angles_df.groupby(["phi_bin", "psi_bin"]).size()

    # Normalizar para obtener probabilidades
    density_table = density_table / density_table.sum()

    return density_table


# Guardar en formato de archivo basado en datos PDB
def save_rama_data(density_table, output_file="pdb_phi_psi_density.cvs"):
    """Guarda la tabla de densidad en un archivo de texto basado en datos de PDB."""
    with open(output_file, "w") as f:
        f.write("# Table name/description: \"PDB-based phi-psi density\"\n")
        f.write("# Number of dimensions: 2\n")
        f.write("# For each dimension, 1 to 2: lower_bound  upper_bound  number_of_bins  wrapping\n")
        f.write("#   x1: -180.0 180.0 180 true\n")
        f.write("#   x2: -180.0 180.0 180 true\n")
        f.write("# List of table coordinates and values. (Value is last number on each line.)\n")

        for (phi, psi), density in density_table.items():
            f.write(f"{phi:.1f} {psi:.1f} {density:.15f}\n")
    print(f"Archivo guardado: {output_file}")


# # Ejecutar el proceso
# density_table = build_phi_psi_density_table("phi_psi_angles_dataset.csv")
# save_rama_data(density_table)


