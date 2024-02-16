import GEOparse
import numpy as np
import pandas as pd
import warnings
import matplotlib
import matplotlib.pyplot as plt
import seaborn

def get(series, class_col='title'):
	gse = GEOparse.get_GEO(geo=series)
	platform_ids = gse.metadata['platform_id']
	samples = gse.phenotype_data.index#['geo_accession']
	classes = gse.phenotype_data['title']
	columns = gse.gsms[samples[0]].table['ID_REF']
	ns = len(samples)
	ng = len(columns)
	ls = []
	for i in range(ns):
		s = gse.gsms[samples[i]].table
		c_drop = s.columns[s.columns != "VALUE"]
		s = s.drop(c_drop, axis=1).transpose()
		ls.append(s)
	df = pd.concat(ls)
	df.index = samples
	df.columns = columns
	df.insert(ng,'CLASS',classes)
	gpls = []
	for pl in platform_ids:
		table = gse.gpls[pl].table
		table.index = table.ID
		gpls.append(table.drop('ID', axis=1))
	if len(gpls)==1:
		gpls = gpls[0]
	else:
		warnings.warn("Sono presenti piattaforme multiple, viene restituita una lista di tabelle gpls")
	return df, gse, gpls
#ESEMPIO DI CHIAMATA VISTO A LEZIONE: df, gse, gpls = geo.get('GSE68849')

def gene_names(cell_ids, gpls, col_name):
	return gpls.loc[cell_ids,col_name].dropna().unique()

def volcano_plot(log2fc, log10pv, thres = (-2,2,2)):
	thx1, thx2, thy = thres
	down = (log2fc < thx1) & (log10pv > thy)
	up = (log2fc > thx2) & (log10pv > thy)
	seaborn.scatterplot(x=log2fc, y=log10pv)
	seaborn.scatterplot(x=log2fc[down], y = log10pv[down])
	seaborn.scatterplot(x=log2fc[up], y = log10pv[up])
	plt.axvline(thx1, color='gray', linestyle='--')
	plt.axvline(thx2, color='gray', linestyle='--')
	plt.axhline(thy, color='gray', linestyle='--')
	plt.title('Volcano Plot')
	plt.xlabel('$log_2 fc$')
	plt.ylabel('$-log_{10} pv$')
	plt.grid()
	plt.show()
	return (down.index[down], up.index[up])

