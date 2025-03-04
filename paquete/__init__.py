# __init__.py

# Exponer funciones desde 'procesamiento_pdb'
from .procesamiento_pdb import process_pdb, get_atom_coords, torsion_angle, classify_residue

# Exponer funciones desde 'graficas'
from .graficas import glycine_ramachandran_plot, proline_ramachandran_plot, general_ramachandran_plot

# Exponer funciones desde 'angles_cal'
from .angles_cal import calculate_angles, calculate_chi1, calculate_chi2

# Exponer funciones desde 'analisis'
from .analisis import calculate_sec_estructure, detect_outliers_by_iqr
