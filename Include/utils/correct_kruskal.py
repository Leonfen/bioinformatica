import scipy
import numpy as np
from pandas import DataFrame
from numpy import ndarray, float32
from typing import Any

def correct_kruskal(G1: DataFrame, G2: DataFrame) ->ndarray[Any, float32]:
    g1_values = G1.values
    g2_values = G2.values
    pvS = []
    for index in range((len(g1_values[0]))):
        try:
            _, pv = scipy.stats.kruskal([g1_values[0][index], g1_values[1][index], g1_values[2][index], g1_values[3][index], g1_values[4][index], g1_values[5][index]], 
                                        [g2_values[0][index], g2_values[1][index], g2_values[2][index], g2_values[3][index], g2_values[4][index], g2_values[5][index]])
            pvS.append(pv)
        except:
            pvS.append(1.0)
    return (np.asarray(pvS))