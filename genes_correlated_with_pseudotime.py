#!/usr/bin/python -i

import pandas as pd
import matplotlib.pyplot as plt
import sys
import numpy as np
import scipy as sc
import sklearn
from sklearn import decomposition
from sklearn import cluster
import seaborn as sns
import math
import IPython
sys.path.insert(0, "/home/hp/CHLA/single-cell-tools/")
from sc_pseudotime import *
from matplotlib.backends.backend_pdf import PdfPages
import os

import argparse

parser = argparse.ArgumentParser(description="runs genes_correlated_with_pseudotime")
parser.add_argument("-e", "--expression-matrix", dest="expr_mat", default="~/single_cell_tools/example_input_files/transcripts.tpm_census_matrix-comma-delimited.csv", help="gene by cell matrix of expression values", metavar="EXPR")
parser.add_argument("-c", "--cell-sets", dest="cell_sets", default="~/single_cell_tools/example_input_files/cell_sets.csv", help="cell sets", metavar="CELL_SETS")
parser.add_argument("-p", "--plot-settings", dest="plot_settings", default="~/single_cell_tools/example_input_files/plot_settings.csv", help="plot settings", metavar="PLOT_SETTINGS")
parser.add_argument("-r", "--corr-method", dest="corr_method", default="spearman", help="method of correlation (spearman or pearson)", metavar="CORR_METHOD")
parser.add_argument("-o", "--outfile", dest="outfile", default="gene_corr_with_ptime", help="a name to give to the output file", metavar="OUTFILE")
parser.add_argument("-pt", "--pseudotime", dest="pseudotime", help="a list of pseudotime values for a set of cells. Can accept multiple values", metavar="PTIME", nargs='+', required=True)


try:
    options = parser.parse_args()
except SystemExit as err: 
	if err.code == 2: 
		parser.print_help()
		sys.exit(0)
 
expression_file = os.path.expanduser(options.expr_mat)
cellset_file    = os.path.expanduser(options.cell_sets)
settings_file   = os.path.expanduser(options.plot_settings)
correlation_method = options.corr_method
out_filename      = options.outfile

pseudotime_files = options.pseudotime

#~ correlation_file = sys.argv[7]


sett = settings(settings_file, cellset_file)
expression_table, annotation = read_expression(expression_file, sett, min_expression = 0.1, min_cells = 5)

## block of code to calculate correlations
pt = map(read_pseudotime_from_file, pseudotime_files)
exit()
correlation = [get_correlation_with_pseudotime(x, expression_table, method=correlation_method) for x in pt]

corr = pd.concat(correlation, axis=1)
corr.columns = pseudotime_files
corr.to_csv("pseudotime_wo_brC_"+correlation_method+"_correlation.csv", sep="\t")
'''
correlation_733 = get_correlation_with_pseudotime(expression_table, pt733, method=correlation_method)
correlation_737 = get_correlation_with_pseudotime(expression_table, pt737, method=correlation_method)
correlation_Ctrl = get_correlation_with_pseudotime(expression_table, ptCtrl, method=correlation_method)
'''
## combine correlations into one DataFrame
'''
corr = pd.DataFrame(0, index=correlation_733.index, columns=["Ctrl","733","737"])
corr["Ctrl"] = correlation_Ctrl
corr["733"]  = correlation_733
corr["737"]  = correlation_737
corr.to_csv("pseudotime_wo_brC/"+correlation_method+"_correlation.csv", sep="\t")
'''
## manually set which correlation method to use
'''
correlation_method = "kendall"
correlation_file = "pseudotime_wo_brC/kendall_correlation.csv"
corr = pd.read_csv(correlation_file, sep="\t", index_col=0)
'''
#~ correlation_method = "spearman"
correlation_file = "pseudotime_wo_brC_spearman_correlation.csv"
corr = pd.read_csv(correlation_file, sep="\t", index_col=0)


## function plots genes of interest (pd.Index) into pdf
def plot_genes_of_interest(genes_of_interest, out_filename, expression_table, annotation, pt):
	pp = PdfPages(out_filename)
	for i,t in enumerate(genes_of_interest):
		fig, ax = plt.subplots(1,3, figsize=(15,5), sharey="row") #define common y axis for set of plots (treatments)
		print i,t
		title = t
		try:
			mg = mygene.MyGeneInfo()
			gene_info = mg.querymany(t, scopes='ensembl.transcript')[0]
			title = t + "  "+ gene_info["symbol"]
			title += "  ("+gene_info["name"]+")"
		except:
			pass
		fig.suptitle(title)
		cntr = 0
		for i in pt:
			plot_gene_with_pseudotime(expression_table, i, t, annotation, ax=ax[0+cntr])
			ax[0+cntr].set_title("733  "+str(cntr)+correlation_method+"=%.2f" % corr.loc[t,pseudotime_files[0+cntr]])
			cntr += 1
		plt.tight_layout()
		plt.subplots_adjust(top=0.85)
		pp.savefig()
		plt.close('all')
	pp.close()

'''
# 733 + 737 - Ctrl
corr["order"] = (corr["733"] + corr["737"] - corr["Ctrl"]).abs()
genes_of_interest = corr.sort_values(by="order", ascending=False).index[:100]
out_filename = "pseudotime_wo_brC/"+correlation_method+"_733+737-Ctrl.pdf"
plot_genes_of_interest(genes_of_interest, out_filename, expression_table, annotation, pt733, pt737, ptCtrl)
'''
# 733
corr["order"] = corr[pseudotime_files[0]].abs()
genes_of_interest = corr.sort_values(by="order", ascending=False).index[:200]
out_filename = "pseudotime_wo_brC_"+correlation_method+"_test.pdf"
plot_genes_of_interest(genes_of_interest, out_filename, expression_table, annotation, pt)
'''
# 737
corr["order"] = corr["737"].abs()
genes_of_interest = corr.sort_values(by="order", ascending=False).index[:200]
out_filename = "pseudotime_wo_brC/"+correlation_method+"_737.pdf"
plot_genes_of_interest(genes_of_interest, out_filename, expression_table, annotation, pt733, pt737, ptCtrl)

# 733 + 737 where (Ctrl has opposite sign) or (abs(Ctrl)< threshold)
threshold = 0.2
corr["order"] = (corr["733"] + corr["737"]).abs()
opposite_sign = corr[ (corr["733"]*corr["Ctrl"] < 0) & (corr["737"]*corr["Ctrl"] < 0)].index
small_abs     = corr[corr["Ctrl"].abs() < threshold].index
genes_of_interest = corr.loc[opposite_sign.union(small_abs)]
genes_of_interest = genes_of_interest.sort_values(by="order", ascending=False).index[:200]
out_filename = "pseudotime_wo_brC/"+correlation_method+"_733+737-where-Ctrl-below-%.2f.pdf" %threshold
plot_genes_of_interest(genes_of_interest, out_filename, expression_table, annotation, pt733, pt737, ptCtrl)

# 733 and 737 above ht and Ctrl below lt or opposite sign
ht = 0.3
lt = 0.2
corr["order"] = (corr["733"] + corr["737"]).abs()
opposite_sign = corr[ (corr["733"]*corr["Ctrl"] < 0) & (corr["737"]*corr["Ctrl"] < 0)].index
small_abs     = corr[corr["Ctrl"].abs() < threshold].index
good_corr_in_knockdown = corr[ (corr["733"].abs() >= ht) & (corr["737"].abs() >= ht) & (corr["733"]*corr["737"] > 0)].index
genes_of_interest = corr.loc[opposite_sign.union(small_abs).intersection(good_corr_in_knockdown)]
genes_of_interest = genes_of_interest.sort_values(by="order", ascending=False).index[:200]
out_filename = "pseudotime_wo_brC/"+correlation_method+"_733-and-737-above-%.2f_Ctrl-below-%.2f.pdf" %(ht,lt)
plot_genes_of_interest(genes_of_interest, out_filename, expression_table, annotation, pt733, pt737, ptCtrl)
'''
#ut = 0.3
#lt = 0.3
#genes_of_interest = corr[(corr["733"]>=ut) & (corr["737"]>=ut) & (corr["Ctrl"]<=lt)]
#genes_of_interest["order"] = genes_of_interest["733"] + genes_of_interest["737"] - genes_of_interest["Ctrl"]
#out_filename = "pseudotime_wo_brC/spearman_733.pdf"
#gene_order = spearman_733["corr"]

if __name__ == "__main__":
	main()
