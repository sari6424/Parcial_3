import pandas as pd

from Parcial_3.paquete import torsion_angle, get_atom_coords



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


