#!/usr/bin/python

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

expression_file = sys.argv[1]
cellset_file    = sys.argv[2]
settings_file   = sys.argv[3]
pseudotime733_file = sys.argv[4]
pseudotime737_file = sys.argv[5]
pseudotimeCtrl_file = sys.argv[6]
correlation_file = sys.argv[7]
correlation_method = sys.argv[8]
out_filename      = sys.argv[9]

sett = settings(settings_file, cellset_file)
expression_table, annotation = read_expression(expression_file, sett, min_expression = 0.1, min_cells = 5)
pt733 = read_pseudotime_from_file(pseudotime733_file)
pt737 = read_pseudotime_from_file(pseudotime737_file)
ptCtrl = read_pseudotime_from_file(pseudotimeCtrl_file)

## block of code to calculate correlations
'''
correlation_method = "spearman"
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
correlation_method = "spearman"
correlation_file = "pseudotime_wo_brC/spearman_correlation.csv"
corr = pd.read_csv(correlation_file, sep="\t", index_col=0)


## function plots genes of interest (pd.Index) into pdf
def plot_genes_of_interest(genes_of_interest, out_filename, expression_table, annotation, pt733, pt737, ptCtrl):
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
		
		plot_gene_with_pseudotime(expression_table, pt733, t, annotation, ax=ax[0])
		ax[0].set_title("733  "+correlation_method+"=%.2f" % corr.loc[t,"733"])
		plot_gene_with_pseudotime(expression_table, pt737, t, annotation, ax=ax[1])
		ax[1].set_title("737  "+correlation_method+"=%.2f" % corr.loc[t,"737"])
		plot_gene_with_pseudotime(expression_table, ptCtrl, t, annotation, ax=ax[2])
		ax[2].set_title("Ctrl  "+correlation_method+"=%.2f" % corr.loc[t,"Ctrl"])
		
		plt.tight_layout()
		plt.subplots_adjust(top=0.85)
		pp.savefig()
		plt.close()
	pp.close()


'''
# 733 + 737 - Ctrl
corr["order"] = (corr["733"] + corr["737"] - corr["Ctrl"]).abs()
genes_of_interest = corr.sort_values(by="order", ascending=False).index[:100]
out_filename = "pseudotime_wo_brC/"+correlation_method+"_733+737-Ctrl.pdf"
plot_genes_of_interest(genes_of_interest, out_filename, expression_table, annotation, pt733, pt737, ptCtrl)
'''
# 733
corr["order"] = corr["733"].abs()
genes_of_interest = corr.sort_values(by="order", ascending=False).index[:200]
out_filename = "pseudotime_wo_brC/"+correlation_method+"_733.pdf"
plot_genes_of_interest(genes_of_interest, out_filename, expression_table, annotation, pt733, pt737, ptCtrl)
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

	
