
import pandas as pd
from paquete.procesamiento_pdb import get_atom_coords
from paquete.procesamiento_pdb import torsion_angle



def calculate_angles(df, get_atom_coords, torsion_angle):
    """
    Calcula los ángulos phi, psi y omega para cada residuo en cada cadena.

    Args:
    - df (pd.DataFrame): DataFrame con los átomos extraídos del PDB.
    - get_atom_coords (function): Función para obtener las coordenadas de un átomo específico.
    - torsion_angle (function): Función para calcular los ángulos de torsión entre 4 átomos.

    Returns:
    - pd.DataFrame: DataFrame con los ángulos phi, psi y omega calculados.
    """
    angles = []

    for chain in df["chain"].unique():
        for residue in df[df["chain"] == chain]["residue num"].unique():
            residue_name = df[(df["chain"] == chain) & (df["residue num"] == residue)]["residue"].values[0]

            # Obtener las coordenadas de los átomos relevantes
            C_ant = get_atom_coords(df, chain, residue - 1, "C")
            N = get_atom_coords(df, chain, residue, "N")
            CA = get_atom_coords(df, chain, residue, "CA")
            C = get_atom_coords(df, chain, residue, "C")
            N_next = get_atom_coords(df, chain, residue + 1, "N")
            C_next = get_atom_coords(df, chain, residue + 1, "CA")

            # Inicializar los ángulos
            phi, psi, omega = float(), float(), float()

            # Calcular phi si todo está bien
            if C_ant is not None and N is not None and CA is not None and C is not None:
                phi = round(torsion_angle(C_ant, N, CA, C), 2)

            # Calcular psi si todo está bien
            if N is not None and CA is not None and C is not None and N_next is not None:
                psi = round(torsion_angle(N, CA, C, N_next), 2)

            # Calcular omega si todo está bien
            if CA is not None and C is not None and N_next is not None and C_next is not None:
                omega = round(torsion_angle(CA, C, N_next, C_next), 2)

            # Guardar los ángulos si no son 0.00
            if phi != 0.00 and psi != 0.00 and omega != 0.00:
                angles.append((chain, residue, residue_name, phi, psi, omega))

    # Crear un DataFrame con los ángulos calculados
    angles_df = pd.DataFrame(angles, columns=["chain", "residue num", "residue", "phi", "psi", "omega"])

    return angles_df

def calculate_chi1(df):
    """
    Calcula los ángulos χ₁ (chi-1) para los residuos de una proteína a partir de las coordenadas atómicas.

    Los ángulos χ₁ son ángulos de torsión entre los átomos N, CA, CB y CG (si están presentes) en residuos de aminoácidos,
    excluyendo glicina y prolina. El cálculo se realiza para cada residuo que tenga los átomos necesarios.

    Args:
    - df (pd.DataFrame): DataFrame que contiene las coordenadas atómicas de los residuos de la proteína. El DataFrame
      debe tener las columnas 'chain', 'residue num', 'residue', 'atom', 'x', 'y', 'z' para representar la estructura.

    Returns:
    - pd.DataFrame: Un DataFrame que contiene los ángulos χ₁ calculados para cada residuo, con las columnas:
      ['chain', 'residue num', 'residue', 'chi1'].
      Si un residuo no tiene los átomos necesarios para calcular χ₁, se excluye del resultado.
    """
    chi1_angles = []

    # Recorre todas las cadenas y residuos
    for chain in df["chain"].unique():
        for residue in df[df["chain"] == chain]["residue num"].unique():
            residue_name = df[(df["chain"] == chain) & (df["residue num"] == residue)]["residue"].values[0]

            if residue_name == "GLY" or residue_name == "PRO":
                continue  # Excluye glicina y prolina

            # Obtén las coordenadas de los átomos necesarios
            coords = get_atom_coords(df, chain, residue, ["N", "CA", "CB"])
            if coords is None:
                continue  # Si falta alguna coordenada, pasa al siguiente residuo

            N, CA, CB = coords

            # Algunas veces los residuos no tienen átomo Cγ, como glicina. Debes manejarlos adecuadamente.
            if "CG" in df[(df["chain"] == chain) & (df["residue num"] == residue)]["atom"].values:
                CG = get_atom_coords(df, chain, residue, "CG")
            else:
                CG = None

            # Si tiene todas las coordenadas necesarias, calcula el ángulo χ₁
            if CG is not None:
                chi1 = torsion_angle(N, CA, CB, CG)
                chi1_angles.append((chain, residue, residue_name, chi1))

    # Crea un DataFrame con los resultados
    chi1_df = pd.DataFrame(chi1_angles, columns=["chain", "residue num", "residue", "chi1"])
    return chi1_df


def calculate_chi2(df):
    """
    Calcula los ángulos χ₂ (chi-2) para los residuos de una proteína a partir de las coordenadas atómicas.

    Los ángulos χ₂ son ángulos de torsión entre los átomos CA, CB, CG (si están presentes) y CD (o átomos equivalentes)
    en residuos de aminoácidos, excluyendo glicina y prolina. El cálculo se realiza para cada residuo que tenga los átomos
    necesarios.

    Args:
    - df (pd.DataFrame): DataFrame que contiene las coordenadas atómicas de los residuos de la proteína. El DataFrame
      debe tener las columnas 'chain', 'residue num', 'residue', 'atom', 'x', 'y', 'z' para representar la estructura.

    Returns:
    - pd.DataFrame: Un DataFrame que contiene los ángulos χ₂ calculados para cada residuo, con las columnas:
      ['chain', 'residue num', 'residue', 'chi2'].
      Si un residuo no tiene los átomos necesarios para calcular χ₂, se excluye del resultado.
    """
    chi2_angles = []

    # Recorre todas las cadenas y residuos
    for chain in df["chain"].unique():
        for residue in df[df["chain"] == chain]["residue num"].unique():
            residue_name = df[(df["chain"] == chain) & (df["residue num"] == residue)]["residue"].values[0]

            if residue_name == "GLY" or residue_name == "PRO":
                continue  # Excluye glicina y prolina

            # Obtén las coordenadas de los átomos necesarios
            coords = get_atom_coords(df, chain, residue, ["CA", "CB"])
            if coords is None:
                continue  # Si falta alguna coordenada, pasa al siguiente residuo

            C_alpha, C_beta = coords

            # Algunas veces los residuos no tienen átomo Cγ o átomos equivalentes.
            C_gamma = None
            if "CG" in df[(df["chain"] == chain) & (df["residue num"] == residue)]["atom"].values:
                C_gamma = get_atom_coords(df, chain, residue, "CG")
            elif "OG" in df[(df["chain"] == chain) & (df["residue num"] == residue)][
                "atom"].values:  # Para serina o treonina
                C_gamma = get_atom_coords(df, chain, residue, "OG")

            # Algunas veces los residuos no tienen átomo Cδ o átomos equivalentes.
            C_delta = None
            if "CD" in df[(df["chain"] == chain) & (df["residue num"] == residue)]["atom"].values:
                C_delta = get_atom_coords(df, chain, residue, "CD")
            elif "OG" in df[(df["chain"] == chain) & (df["residue num"] == residue)][
                "atom"].values:  # Para serina o treonina
                C_delta = get_atom_coords(df, chain, residue, "OG")

            # Si no se pueden obtener las coordenadas de todos los átomos, salta este residuo
            if C_alpha is None or C_beta is None or C_gamma is None or C_delta is None:
                continue  # Salta este residuo y pasa al siguiente

            # Si tiene todas las coordenadas necesarias, calcula el ángulo χ₂
            try:
                chi2 = torsion_angle(C_alpha, C_beta, C_gamma, C_delta)
                chi2_angles.append((chain, residue, residue_name, chi2))
            except Exception:
                continue  # En caso de que no se pueda calcular el ángulo, también se salta este residuo

    # Crea un DataFrame con los resultados
    chi2_df = pd.DataFrame(chi2_angles, columns=["chain", "residue num", "residue", "chi2"])
    return chi2_df

