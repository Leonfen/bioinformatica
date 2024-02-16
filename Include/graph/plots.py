from pandas import DataFrame
import numpy as np
from utils import geo

def plotting(sample1: int, sample2: int, pv, df: DataFrame, treshold: tuple = (-2, 2, 2)):
    log10pv = -1 * np.log10(pv)  # Asse y
    means = df.groupby('CLASS').mean()  # Medie fra le classi, fra la classe influenza e virus, rispetto ai gruppi in
    # cui ha diviso i campioni
    dm = means.iloc[sample1] / means.iloc[sample2]  # rapporto fra le medie, mean.iloc[0] sarebbe rapporto tra il primo
    # e secondo campione
    log2dm = np.log2(dm)  # Asse x, logaritmo in base due fra i rapporti delle medie
    return geo.volcano_plot(log2dm, log10pv, treshold)

