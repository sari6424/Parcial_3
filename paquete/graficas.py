import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def glycine_ramachandran_plot(angles_df):
    """
    Genera un gráfico de Ramachandran para residuos de glicina.
    Filtra los residuos de glicina antes de generar el gráfico.

    Args:
    - angles_df (pd.DataFrame): DataFrame con las columnas 'phi', 'psi', y 'residue'.

    Returns:
    - plt.figure: El objeto figura del gráfico.
    """
    # Filtra los residuos específicos de glicina
    glycine_angles_df = angles_df.loc[angles_df['residue'] == 'GLY']

    plt.figure(figsize=(8, 6))

    # Crear el gráfico para glicina
    sns.scatterplot(x="phi", y="psi", data=glycine_angles_df, s=15, edgecolor="green", color="lightgreen")
    sns.kdeplot(x="phi", y="psi", data=glycine_angles_df, levels=3, color="darkgreen", fill=True, alpha=0.3,
                bw_adjust=0.25)

    # Configuración de la gráfica
    plt.xlabel("Phi (φ)")
    plt.ylabel("Psi (ψ)")
    plt.xlim(-180, 180)
    plt.ylim(-180, 180)
    plt.title("Ramachandran plot for glycine residues")
    plt.grid(True)

    return plt


def proline_ramachandran_plot(angles_df):
    """
    Genera un gráfico de Ramachandran para residuos de prolina.
    Filtra los residuos de prolina antes de generar el gráfico.

    Args:
    - angles_df (pd.DataFrame): DataFrame con las columnas 'phi', 'psi', y 'residue'.

    Returns:
    - plt.figure: El objeto figura del gráfico.
    """
    # Filtra los residuos específicos de prolina
    proline_angles_df = angles_df.loc[angles_df['residue'] == 'PRO']

    plt.figure(figsize=(8, 6))

    # Crear el gráfico para prolina
    sns.scatterplot(x="phi", y="psi", data=proline_angles_df, s=15, edgecolor="red", color="tomato")
    sns.kdeplot(x="phi", y="psi", data=proline_angles_df, levels=3, color="darkred", fill=True, alpha=0.3,
                bw_adjust=0.3)

    # Configuración de la gráfica
    plt.xlabel("Phi (φ)")
    plt.ylabel("Psi (ψ)")
    plt.xlim(-180, 180)
    plt.ylim(-180, 180)
    plt.title("Ramachandran plot for proline residues")
    plt.grid(True)

    return plt


def general_ramachandran_plot(angles_df):
    """
    Genera un gráfico de Ramachandran para todos los residuos.

    Args:
    - angles_df (pd.DataFrame): DataFrame con las columnas 'phi' y 'psi'.

    Returns:
    - plt.figure: El objeto figura del gráfico.
    """
    plt.figure(figsize=(8, 6))

    # Crear el gráfico general
    sns.scatterplot(x="phi", y="psi", data=angles_df, s=15, edgecolor="blue", color="lightblue")
    sns.kdeplot(x="phi", y="psi", data=angles_df, levels=4, color="darkblue", fill=True, alpha=0.3)

    # Configuración de la gráfica
    plt.xlabel("Phi (φ)")
    plt.ylabel("Psi (ψ)")
    plt.xlim(-180, 180)
    plt.ylim(-180, 180)
    plt.title("Ramachandran plot")
    plt.grid(True)

    return plt