# Ramachandran package

Este paquete en Python permite procesar archivos PDB para calcular los ángulos de torsión phi, psi y omega de una estructura proteica, analizar su estructura secundaria y detectar outliers. Además, proporciona herramientas para visualizar los ángulos en gráficos de Ramachandran.

## Características
- Procesamiento de archivos PDB para extraer coordenadas atómicas.
- Cálculo de los ángulos de torsión phi, psi y omega.
- Análisis de estructura secundaria basado en helices y láminas beta.
- Detección de valores atípicos en los ángulos de torsión.
- Generación de gráficos de Ramachandran para glicina, prolina y otros residuos.

## Instalación
Este paquete requiere Python 3.x y las siguientes dependencias:

```bash
pip install numpy pandas matplotlib seaborn
```

Además, asegúrate de tener BLAST+ instalado si deseas verificar la especificidad de las secuencias.

## Uso

### Procesamiento de un archivo PDB

```python
from paquete.procesamiento_pdb import process_pdb

data = process_pdb("estructura.pdb")
print(data.head())
```

### Cálculo de ángulos de torsión

```python
from paquete.angles_cal import calculate_angles

angles_df = calculate_angles(data)
print(angles_df.head())
```

### Análisis de estructura secundaria

```python
from paquete.analisis import calculate_sec_estructure

helices, laminas, total = calculate_sec_estructure("estructura.pdb")
print(f"Hélices: {helices}%, Láminas: {laminas}%, Total residuos: {total}")
```

### Detección de valores atípicos en los ángulos

```python
from paquete.analisis import detect_outliers_by_iqr

angles_outliers = detect_outliers_by_iqr(angles_df)
print(angles_outliers.head())
```

### Generación de gráficos de Ramachandran

```python
from paquete.graficas import general_ramachandran_plot

general_ramachandran_plot(angles_df).show()
```


## Licencia
Este proyecto está licenciado bajo la MIT License.

