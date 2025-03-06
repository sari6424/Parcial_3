# __init__.py

# Exponer funciones desde 'procesamiento_pdb'
from .procesamiento_pdb import process_pdb, get_atom_coords, torsion_angle, classify_residue, save_angles

# Exponer funciones desde 'graficas'
from .graficas import glycine_ramachandran_plot, proline_ramachandran_plot, general_ramachandran_plot

# Exponer funciones desde 'angles_cal'
from .angles_cal import calculate_angles

# Exponer funciones desde 'analisis'
from .analisis import calculate_sec_estructure, detect_outliers_by_iqr

# Exponer funciones desde 'density_data'
from .density_data import build_phi_psi_density_table, save_rama_data, load_phi_psi_data

# Exponer funciones desde 'make_data_base'
from .make_data_base import create_pdb_folder, download_pdb, download_pdbs_from_txt, process_multiple_pdbs
