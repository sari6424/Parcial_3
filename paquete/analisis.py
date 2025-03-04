def calculate_sec_estructure(pdb_file):
    """
    Calcula el porcentaje de residuos en hélices alfa y láminas beta a partir de un archivo PDB.

    Parámetros:
    - pdb_file (str): Ruta del archivo PDB.

    Retorna:
    - porcentaje_helice (float): Porcentaje de residuos en hélices alfa.
    - porcentaje_lamina (float): Porcentaje de residuos en láminas beta.
    """
    helix = 0
    sheet = 0
    tot_residu = 0

    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith("ATOM") and line[13:15] == "CA":  # Contamos Cα (residuos totales)
                tot_residu += 1
            elif line.startswith("HELIX"):  # Hélice alfa
                helix += int(line[71:76].strip())  # Número de residuos en la hélice
            elif line.startswith("SHEET"):  # Lámina beta
                try:
                    start = int(line[23:26].strip())  # Residuo inicial de la lámina
                    end = int(line[34:37].strip())  # Residuo final de la lámina
                    sheet += (end - start + 1)  # Cantidad de residuos en la lámina
                except ValueError:
                    continue  # Evita errores si hay problemas con los valores

    # Cálculo de porcentajes con prevención de división por cero
    percent_helix = (helix / tot_residu) * 100 if tot_residu else 0
    percent_sheet = (sheet / tot_residu) * 100 if tot_residu else 0

    return percent_helix, percent_sheet, tot_residu

def detect_outliers_by_iqr(angles_df):
    """
    Detecta outliers en las columnas 'phi' y 'psi' del DataFrame usando el método IQR
    La función buscará las columnas con nombre 'phi' y 'psi' en el DataFrame y marcará los outliers.

    Args:
    - angles_df (pd.DataFrame): DataFrame que puede contener varias columnas, pero debe tener las columnas 'phi' y 'psi'.

    Returns:
    - pd.DataFrame: El DataFrame con nuevas columnas 'phi_is_outlier' y 'psi_is_outlier' que indican si el valor es un outlier.
      Además, se agrega una columna 'is_outlier' que marca si hay outliers en cualquiera de las dos columnas.
    """

    # Verificar que las columnas 'phi' y 'psi' existen en el DataFrame
    if 'phi' not in angles_df.columns or 'psi' not in angles_df.columns:
        raise ValueError("El DataFrame debe contener las columnas 'phi' y 'psi'")

    for angle in ["phi", "psi"]:
        # Calcular los cuartiles y el rango intercuartílico
        Q1 = angles_df[angle].quantile(0.25)
        Q3 = angles_df[angle].quantile(0.75)
        IQR = Q3 - Q1
        lower_bound = Q1 - 1.5 * IQR
        upper_bound = Q3 + 1.5 * IQR

        # Marcar los outliers
        angles_df[f"{angle}_is_outlier"] = ~angles_df[angle].between(lower_bound, upper_bound)

    # Crear una columna para indicar si un residuo tiene algún outlier en phi o psi
    angles_df["is_outlier"] = angles_df["phi_is_outlier"] | angles_df["psi_is_outlier"]

    return angles_df


