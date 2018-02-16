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
import os
import ipdb 

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

pseudotime_files = sorted(options.pseudotime)

#~ correlation_file = sys.argv[7]


sett = settings(settings_file, cellset_file)
expression_table, annotation = read_expression(expression_file, sett, min_expression = 0.1, min_cells = 5)

#~ correlation_method = "spearman"



#~ correlation_file = "pseudotime_wo_brC_spearman_correlation.csv"
#~ corr = pd.read_csv(correlation_file, sep="\t", index_col=0)
pt = map(read_pseudotime_from_file, pseudotime_files)
ptime_titles = [i.replace(".csv", "").rsplit("/")[-1] for i in pseudotime_files]
pt = dict(zip(ptime_titles, pt))

correlation_file = "pseudotime_wo_brC_spearman_correlation.csv"
corr = pd.read_csv(correlation_file, sep="\t", index_col=0)

def set_pts(pseudotime_files):
	pt = map(read_pseudotime_from_file, pseudotime_files)
	correlation = [get_correlation_with_pseudotime(x, expression_table, method=correlation_method) for x in pt]
	corr = pd.concat(correlation, axis=1)
	corr.columns = pt.keys()
	ptime_titles = [i.replace(".csv", "").rsplit("/")[-1] for i in pseudotime_files]
	pt = dict(zip(ptime_titles, pt))
	corr.to_csv("pseudotime_wo_brC_"+correlation_method+"_correlation.csv", sep="\t")
	return corr, pt


if sorted(pt.keys()) == sorted(corr.columns):
	pass
else:
	corr, pt = set_pts(pseudotime_files)
	corr.columns = sorted(pt.keys())
	corr.to_csv("pseudotime_wo_brC_"+correlation_method+"_correlation.csv", sep="\t")

#~ IPython.embed()

user_ptimes = ' '.join(pt.keys())

## function plots genes of interest (pd.Index) into pdf
def plot_genes_of_interest(genes_of_interest, out_filename, expression_table, annotation, pt_dict):
	out_genes = out_filename.replace(".pdf", "_genes.csv")
	og = open(out_genes, 'w')
	pp = PdfPages(out_filename)
	for i,t in enumerate(genes_of_interest):
		fig, ax = plt.subplots(1,len(pt_dict), figsize=(15,5), sharey="row") #define common y axis for set of plots (treatments)
		print i,t
		title = t
		try:
			mg = mygene.MyGeneInfo()
			gene_info = mg.querymany(t, scopes='ensembl.transcript')[0]
			title = t + "  "+ gene_info["symbol"]
			title += "  ("+gene_info["name"]+")"
			og.write(title+"\n")
			#~ print corr.loc[t,pt_dict.keys()[0]]
		except:
			pass
		fig.suptitle(title)
		cntr = 0
		for k,v in pt_dict.iteritems():
			plot_gene_with_pseudotime(expression_table, v, t, annotation, ax=ax[0+cntr])
			ax[0+cntr].set_title(k+correlation_method+"=%.2f" % corr.loc[t,k])
			cntr += 1
			#~ plt.show()
		plt.tight_layout()
		plt.subplots_adjust(top=0.85)
		pp.savefig()
	pp.close()
	og.close()
	
		#~ plt.tight_layout()
		#~ plt.subplots_adjust(top=0.85)
		#~ pp.savefig()
		#~ plt.close('all')
	#~ pp.close()

'''
# 733 + 737 - Ctrl
corr["order"] = (corr["733"] + corr["737"] - corr["Ctrl"]).abs()
genes_of_interest = corr.sort_values(by="order", ascending=False).index[:100]
out_filename = "pseudotime_wo_brC/"+correlation_method+"_733+737-Ctrl.pdf"
plot_genes_of_interest(genes_of_interest, out_filename, expression_table, annotation, pt733, pt737, ptCtrl)
'''

# main loop / choosing action
while True:
	question = """Choose from following:
	[C]	Create Correlation file from New Pseudotime Files
	[D]	Plot User-Supplied Genes
	[T]	Plot Top N Genes with highest correlation
	[X]	Exit
	"""
	action = raw_input(question).upper()
	if(action == "X"):
		break
		#~ exit()
		#~ FACS_0407_2017_SHL_input_files/DEGS_day_12.csv
	elif(action == "C"):
		ptime_paths = raw_input("provide path(s) to new pseudotime files ").split(",")
		corr_out = raw_input("provide filename for new correlation files ")
		## block of code to calculate correlations
		pt = map(read_pseudotime_from_file, ptime_paths)
		corr = [get_correlation_with_pseudotime(x, expression_table, method=correlation_method) for x in pt]
		corr.to_csv(corr_out+"_"+correlation_method+"_correlation.csv", sep="\t")
		ptime_titles = [i.replace(".csv", "").rsplit("/")[-1] for i in ptime_paths]
		ptimes = dict(zip(ptime_titles, ptime_paths))
		user_ptimes = ' '.join(ptime_titles)
		corr.columns = ptime_titles
		
	elif(action == "D"):
		DEG_path = raw_input("provide path to differentially expressed genes ")
		ptime = raw_input("Which pseudotime would you like correlate with? ("+user_ptimes+ ") ")
		DEGS = pd.read_csv(DEG_path, index_col=0)
		corr["order"] = corr[ptime].abs()
		DEGS = corr[corr.index.isin(DEGS.index)].index
		out_filename = "pseudotime_wo_brC/"+correlation_method+"_"+ptime+"_DEGS.pdf"
		plot_genes_of_interest(DEGS, out_filename, expression_table, annotation, pt)
	elif(action == "T"):
		top_n = int(raw_input("How many genes would you like to plot? "))
		ptime = raw_input("Which pseudotime would you like correlate with? ("+user_ptimes+ ") ")
		# 733
		corr["order"] = corr[ptime].abs()
		genes_of_interest = corr.sort_values(by="order", ascending=False).index[:top_n]
		out_filename = "pseudotime_wo_brC/"+correlation_method+"_"+ptime+".pdf"
		plot_genes_of_interest(genes_of_interest, out_filename, expression_table, annotation, pt)
		print genes_of_interest
	elif(action == "I"):
		IPython.embed()

'''
#DEGS
DEGS = pd.read_csv("FACS_0407_2017_SHL_input_files/DEGS_day_12.csv", index_col=0)
corr["order"] = corr[ptime_titles[0]].abs()
DEGS = corr[corr.index.isin(DEGS.index)].index
out_filename = "pseudotime_wo_brC/"+correlation_method+"_733.pdf"
plot_genes_of_interest(DEGS, out_filename, expression_table, annotation, pt)

# 733
corr["order"] = corr[ptime_titles[0]].abs()
genes_of_interest = corr.sort_values(by="order", ascending=False).index[:200]
out_filename = "pseudotime_wo_brC/"+correlation_method+"_733_20180214.pdf"
plot_genes_of_interest(genes_of_interest, out_filename, expression_table, annotation, pt)

# 737
corr["order"] = corr[ptime_titles[1]].abs()
genes_of_interest = corr.sort_values(by="order", ascending=False).index[:200]
out_filename = "pseudotime_wo_brC/"+correlation_method+"_737.pdf"
plot_genes_of_interest(genes_of_interest, out_filename, expression_table, annotation, pt)

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

#~ if __name__ == "__main__":
	#~ main()
