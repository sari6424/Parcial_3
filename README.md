# Ramachandran Package

Este paquete en Python permite procesar archivos PDB para calcular los ángulos de torsión **phi (ϕ), psi (ψ) y omega (ω)** de una estructura proteica, analizar su estructura secundaria, detectar valores atípicos y generar diagramas de Ramachandran. Además, proporciona herramientas para calcular la densidad de ocupación de ángulos Ramachandran a partir de datos PDB.

## Características

- Procesamiento de archivos **PDB** para extraer coordenadas atómicas.
- Cálculo de los ángulos de torsión **phi (ϕ), psi (ψ) y omega (ω)**.
- Análisis de estructura secundaria basado en **hélices alfa y láminas beta**.
- Detección de valores atípicos en los ángulos de torsión.
- Generación de gráficos de Ramachandran.
- Cálculo de la densidad de ocupación de ángulos Ramachandran con datos PDB.
- Descarga automática de archivos PDB.

---

## Instalación

Este paquete requiere **Python 3.x** y las siguientes dependencias:

```bash
pip install numpy pandas matplotlib seaborn requests os
```

---

## Uso

### 1. Procesamiento de un archivo PDB

```python
from paquete.procesamiento_pdb import process_pdb

# Procesar un único archivo PDB
data = process_pdb("estructura.pdb")
print(data)
```

---

### 2. Procesamiento de múltiples archivos PDB

```python
from paquete.make_data_base import process_multiple_pdbs

# Procesar una carpeta con múltiples archivos PDB
angles_data = process_multiple_pdbs("pdb_files/")
print(angles_data.head())
```

---

### 3. Cálculo de ángulos de torsión

```python
from paquete.angles_cal import calculate_angles

# Calcular ángulos phi, psi y omega de una estructura PDB
angles_df = calculate_angles(data)
print(angles_df)
```

---

### 4. Análisis de estructura secundaria

```python
from paquete.analisis import calculate_sec_estructure

# Analizar la estructura secundaria de una proteína
helices, laminas, total = calculate_sec_estructure("estructura.pdb")
print(f"Hélices: {helices}%, Láminas: {laminas}%, Total residuos: {total}")
```

---

### 5. Detección de valores atípicos en los ángulos

```python
from paquete.analisis import detect_outliers_by_iqr

# Detectar valores atípicos en los ángulos de torsión
angles_outliers = detect_outliers_by_iqr(angles_df)
print(angles_outliers)
```

---

### 6. Generación de gráficos de Ramachandran

```python
from paquete.graficas import general_ramachandran_plot

# Generar el gráfico de Ramachandran
general_ramachandran_plot(angles_df).show()
```

---

### 7. Cálculo de la densidad de ocupación de ángulos Ramachandran

```python
from paquete.density_data import build_phi_psi_density_table, save_rama_data

# Construir tabla de densidad
density_table = build_phi_psi_density_table("phi_psi_angles_dataset.csv")

# Guardar la tabla en un archivo
save_rama_data(density_table, "pdb_phi_psi_density.data")
```

---

### 8. Descarga automática de archivos PDB

```python
from paquete.make_data_base import download_pdb, download_pdbs_from_txt

# Descargar un solo PDB
download_pdb("1CRN")

# Descargar múltiples PDBs desde un archivo TXT
download_pdbs_from_txt("listapdb.txt")
```


---

## Licencia

Este proyecto está bajo la licencia MIT.



