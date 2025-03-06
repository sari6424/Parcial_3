import numpy as np
import pandas as pd


def torsion_angle(a, b, c, d):
    """
    Calcula el ángulo de torsión entre cuatro átomos a, b, c, d.

    Args:
    - a, b, c, d (array-like): Coordenadas (x, y, z) de los átomos.

    Returns:
    - float: El ángulo de torsión en grados.
    """
    # Calcula vectores entre los puntos
    u1 = np.array(b) - np.array(a)  # Vector entre a y b
    u2 = np.array(c) - np.array(b)  # Vector entre b y c
    u3 = np.array(d) - np.array(c)  # Vector entre c y d

    # Calcula productos cruzados y normaliza para el ángulo
    e = np.linalg.norm(u2) * u1
    f = np.cross(u2, u3)  # Producto cruzado de los vectores b-c y c-d
    g = np.cross(u1, u2)  # Producto cruzado de los vectores a-b y b-c

    # Calcula los componentes del ángulo
    y = np.dot(e, f)  # Producto punto entre e y f
    x = np.dot(g, f)  # Producto punto entre g y f

    # Calcula y retorna el ángulo en grados
    ang = np.arctan2(y, x)  # Función atan2 para calcular el ángulo
    return np.degrees(ang)  # Convertir a grados

def get_atom_coords(df, chain, residue, atom):
    """
    Obtiene las coordenadas (x, y, z) de un átomo específico en un DataFrame.

    Args:
    - df (pd.DataFrame): DataFrame que contiene los datos de los átomos.
    - chain (str): La cadena del residuo.
    - residue (int): El número de residuo.
    - atom (str): El nombre del átomo (por ejemplo, "N", "CA", "C").

    Returns:
    - list or None: Las coordenadas (x, y, z) del átomo o None si no se encuentra.
    """
    try:
        # Filtra el DataFrame para obtener la fila del átomo solicitado
        atom_row = df[(df["chain"] == chain) & (df["residue num"] == residue) & (df["atom"] == atom)]
        return atom_row[['x', 'y', 'z']].values[0]  # Devuelve las coordenadas del átomo
    except IndexError:
        return None  # Si no se encuentra el átomo, retorna None



def classify_residue(residue):
    """
    Clasifica un residuo de aminoácido en una de las categorías: polar, polar cargado positivo,
    polar cargado negativo, no polar alifático, o no polar aromático.

    Args:
    - residue (str): El nombre del residuo de aminoácido (por ejemplo, "ASN", "CYS").

    Returns:
    - str: La categoría del residuo (por ejemplo, "Polar", "Polar positive charged").
    """
    # Definir listas de categorías de residuos
    polar = ["ASN", "CYS", "GLN", "SER", "THR"]
    polar_pos = ["ARG", "HIS", "LYS"]
    polar_neg = ["ASP", "GLU"]
    nonpolar_ali = ["ALA", "ILE", "GLY", "LEU", "MET", "PRO", "VAL"]
    nonpolar_aro = ["PHE", "TYR", "TRP"]

    # Clasificación basada en la lista correspondiente
    if residue in polar:
        return "Polar"
    elif residue in polar_pos:
        return "Polar positive charged"
    elif residue in polar_neg:
        return "Polar negative charged"
    elif residue in nonpolar_ali:
        return "Non polar aliphatic"
    elif residue in nonpolar_aro:
        return "Non polar aromatic"
    else:
        return "Unknown residue"  # Si no se encuentra el residuo en ninguna categoría


def process_pdb(pdb_file):
    """
    Procesa un archivo PDB y extrae los átomos que corresponden a aminoácidos.

    Args:
    - pdb_file (str): Ruta al archivo PDB.

    Returns:
    - pd.DataFrame: DataFrame con la información de los átomos extraídos.
    """
    data = []
    # Lista de aminoácidos válidos
    aminoacids = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
                  "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]

    # Abre el archivo PDB y procesa las líneas
    with open(pdb_file) as file:
        for line in file:
            if line.startswith("ATOM"):
                atom = line[12:16].strip()
                residue = line[17:20].strip()
                chain = line[21].strip()
                res_num = int(line[22:26].strip())
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())

                # Incluye solo los átomos de aminoácidos válidos
                if residue in aminoacids:
                    data.append([atom, residue, chain, res_num, x, y, z])

    # Crear un DataFrame con los datos extraídos
    titles = ["atom", "residue", "chain", "residue num", "x", "y", "z"]
    return pd.DataFrame(data, columns=titles)


def save_angles(angles_df, output_file):
    """
    Guarda el DataFrame con los ángulos phi, psi y omega en un archivo TSV.

    Args:
    - angles_df (pd.DataFrame): DataFrame con los ángulos calculados.
    - output_file (str): Nombre del archivo de salida (por defecto "angles_data.tsv").

    Returns:
    - None
    """
    angles_df.to_csv(output_file, sep="\t", index=False)
    print(f"Archivo guardado exitosamente en: {output_file}")
