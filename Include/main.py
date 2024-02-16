import scipy as sc
from utils import geo
from graph import plots
from graph.tests import kruskal, t_test, oneway, utest
from utils.output import all_gene_names

def main():
    df, gse, gpls = geo.get('GSE66597')
    ns, ng = df.shape  # take the number of samples and genes
    df.iloc[[0, 1, 2, 3, 4, 5, 6, 7, 8], ng-1] = 'Control'  # Gene's expressed levels in after 0h
    df.iloc[[9, 10, 11, 12, 13, 14, 15, 16, 17], ng-1] = 'Infected'  # Expressed levels after 24h, exposed on LPS and LbTR enzymes
    # Take all the samples
    G1 = df.loc[df.CLASS== 'Control', df.columns != "CLASS"] # After 0h 
    G2 = df.loc[df.CLASS == 'Infected', df.columns != "CLASS"] # After 24h
    
    t_test_pv, t_test_gene_names = t_test(G1, G2, gpls, 0.000003) 
    oneway_pv, oneway_gene_names= oneway(G1, G2, gpls, 0.0000028)
    kurskal_pv, kruskal_gene_names = kruskal(G1, G2, gpls, 0.00394) 
    utest_pv, utest_gene_names = utest(G1, G2, gpls, 0.00041229480206169137)
    
    t_test_pv_manual, t_test_gene_names_manual = t_test(G1, G2, gpls, 0.00000002, True)
    oneway_pv_manual, oneway_gene_names_manual = oneway(G1, G2, gpls, 0.000000021, True)
    kurskal_pv_manual, kruskal_gene_names_manual = kruskal(G1, G2, gpls, 0.0038849066221194816, True)
    utest_pv_manual, utest_gene_names_manual = utest(G1, G2, gpls, 0.0004099601088531576, True)


    # Vulcano plot - normal
    down, up = plots.plotting(0, 1, t_test_pv, df, (-0.65, 0.65, 3.0))
    # Overexpressed and underexpressed genes
    gene_names_down = geo.gene_names(down, gpls, 'GENE_SYMBOL') 
    gene_names_up = geo.gene_names(up, gpls, 'GENE_SYMBOL')
    
    all_gene_names('Over Expressed Genes', gene_names_up)
    all_gene_names('Under Expressed Genes', gene_names_down)

    # Vulcano plot - manual
    down, up = plots.plotting(0, 1, t_test_pv, df, (-0.82, 0.82, 3.5))
    # Overexpressed and underexpressed genes
    gene_names_down = geo.gene_names(down, gpls, 'GENE_SYMBOL') 
    gene_names_up = geo.gene_names(up, gpls, 'GENE_SYMBOL')
    
    all_gene_names('Over Expressed Genes', gene_names_up, True)
    all_gene_names('Under Expressed Genes', gene_names_down, True)

main()