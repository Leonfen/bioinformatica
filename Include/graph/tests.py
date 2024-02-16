import scipy
from pandas import DataFrame
from utils import geo
from utils.output import all_gene_names
from numpy import ndarray, float32
from typing import Any, Tuple, List
from utils.correct_kruskal import correct_kruskal

def t_test(G1: DataFrame, G2: DataFrame, gpls: DataFrame, alpha: float = 0.008, manual: bool = False, correct_pvalue: bool = False) -> Tuple[ndarray[Any, float32], List[str]]:
    pv: ndarray[Any, float32] = scipy.stats.ttest_ind(G1, G2).pvalue # Prendi il pvalue
    col_Test: Any = G1.columns[pv < alpha] # Prendi tutti i geni con un pvalue minore di alpha
    gene_names: List[str] = geo.gene_names(col_Test, gpls, "GENE_SYMBOL") # Prendi tutti i nomi dei geni significativi
    test = 'T-test genes' if not manual else 'T-test manual genes'
    if correct_pvalue: pv = scipy.stats.false_discovery_control(pv) # Correct p-value
    all_gene_names(test, gene_names, manual)
    return (pv, gene_names)

def utest(G1: DataFrame, G2: DataFrame, gpls: DataFrame, alpha: float = 0.008, manual: bool = False) -> Tuple[ndarray[Any, float32], List[str]]:
    pv: ndarray[Any, float32] = scipy.stats.mannwhitneyu(G1, G2).pvalue
    col_utest = G1.columns[pv < alpha]
    gene_names: List[str] = geo.gene_names(col_utest, gpls, "GENE_SYMBOL") # Get names of all significant genes
    test = 'Utest genes' if not manual else 'U-test manual genes'
    all_gene_names(test, gene_names, manual)
    return (pv, gene_names)
    

def oneway(G1: DataFrame, G2: DataFrame, gpls: DataFrame, alpha: float = 0.008, manual: bool = False) -> Tuple[ndarray[Any, float32], List[str]]:
    pv: ndarray[Any, float32] = scipy.stats.f_oneway(G1, G2).pvalue
    col_anova = G1.columns[pv < alpha]
    gene_names: List[str] = geo.gene_names(col_anova, gpls, "GENE_SYMBOL") # Get names of all significant genes
    test = 'Oneway genes' if not manual else 'Oneway manual genes'
    all_gene_names(test, gene_names, manual)
    return (pv, gene_names)


def kruskal(G1: DataFrame, G2: DataFrame, gpls: DataFrame, alpha: float = 0.008, manual: bool = False) -> Tuple[ndarray[Any, float32], List[str]]:
    pv: ndarray[Any, float32] = correct_kruskal(G1, G2).pvalue
    col_kruskal = G1.columns[pv < alpha]
    gene_names: List[str] = geo.gene_names(col_kruskal, gpls, "GENE_SYMBOL") # Get names of all significant genes
    test = 'Kruskal genes' if not manual else 'Kruskal manual genes'
    all_gene_names(test, gene_names, manual)
    return (pv, gene_names)
