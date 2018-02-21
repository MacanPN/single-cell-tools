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
parser.add_argument("-pt", "--pseudotime", dest="pseudotime", help="experimental cells. a list of pseudotime values. Can accept multiple values", metavar="PTIME", nargs='+', required=True)
parser.add_argument("-cpt", "--control-pseudotime", dest="ctrl_pseudotime", help="control cells. a list of pseudotime values. Can accept multiple values", metavar="PTIME", nargs='+', required=True)


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
output_dir      = options.outfile+"/"

pseudotime_files = sorted(options.pseudotime)
ctrl_pseudotime_files = sorted(options.ctrl_pseudotime) 

sett = settings(settings_file, cellset_file)
expression_table, annotation = read_expression(expression_file, sett, min_expression = 0.1, min_cells = 5)

# read in pseudotime files
pt = map(read_pseudotime_from_file, pseudotime_files)
ptime_titles = [i.replace(".csv", "").rsplit("/")[-1] for i in pseudotime_files]
pt = dict(zip(ptime_titles, pt))

cpt = map(read_pseudotime_from_file, ctrl_pseudotime_files)
cptime_titles = [i.replace(".csv", "").rsplit("/")[-1] for i in ctrl_pseudotime_files]
cpt = dict(zip(cptime_titles, cpt))

correlation_file = "_".join(ptime_titles)+"_"+correlation_method+"_correlation.csv"
ctrl_correlation_file = "_".join(cptime_titles)+"_"+correlation_method+"_correlation.csv"

# read correlation files from similarly named files
if os.path.exists(correlation_file):
	corr = pd.read_csv(correlation_file, sep="\t", index_col=0)
if os.path.exists(ctrl_correlation_file):
	ctrl_corr = pd.read_csv(ctrl_correlation_file, sep="\t", index_col=0)
 
def set_pts(pseudotime_files, cell_set_flag):
	pt = map(read_pseudotime_from_file, pseudotime_files)
	correlation = [get_correlation_with_pseudotime(x, expression_table, annotation, cell_set_flag, method=correlation_method) for x in pt]
	corr = pd.concat(correlation, axis=1)
	ptime_titles = [i.replace(".csv", "").rsplit("/")[-1] for i in pseudotime_files]
	pt = dict(zip(ptime_titles, pt))
	return corr, pt

corr_columns = []
for i in sorted(pt.keys()):
	corr_columns += [i+"_exp_corr"]
	corr_columns += [i+"_ctrl_corr"]

# check if experimental correlation files have already been read in
try:
	corr
except:
	corr, pt = set_pts(pseudotime_files, cell_set_flag="exp")
	corr_columns = []
	for i in sorted(pt.keys()):
		corr_columns += [i+"_exp_corr"]
		corr_columns += [i+"_ctrl_corr"]
	corr.columns = corr_columns
	corr.to_csv(correlation_file, sep="\t")
else:
	if corr_columns == list(corr.columns):
		pass
	else:
		print "column names do not match!"

# check if control correlation files have already been read in
try:
	ctrl_corr
except:
	ctrl_corr, cpt = set_pts(ctrl_pseudotime_files, cell_set_flag="ctrl")
	ctrl_corr.columns = sorted(cpt.keys())

	ctrl_corr.to_csv(ctrl_correlation_file, sep="\t")
else:
	if sorted(cpt.keys()) == list(ctrl_corr.columns):
		pass
	else:
		print "column names do not match!"

user_ptimes = ' '.join(pt.keys())
ctrl_user_ptimes = ' '.join(cpt.keys())

## function plots genes of interest (pd.Index) into pdf
def plot_genes_of_interest(genes_of_interest, out_filename, expression_table, annotation, pt, ctrl_pseudotime=None):
	#~ IPython.embed()
	plot_id = ["RBKD", "Ctrl_wo_RBKD", "Ctrl_alone"]
	out_genes = out_filename.replace(".pdf", "_genes.csv")
	og = open(out_genes, 'w')
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
			og.write(title+"\n")
			#~ print corr.loc[t,pt_dict.keys()[0]]
		except:
			pass
		fig.suptitle(title)
		cntr = 0
		#~ for k,v in pt_dict.iteritems():
		while cntr < 3:
			plot_gene_with_pseudotime(expression_table, pt, t, annotation, ax=ax[0+cntr], plot_id=plot_id[cntr], ctrl_pseudotime=ctrl_pseudotime)
			#~ ax[0+cntr].set_title(ptime+correlation_method+"=%.2f" % corr.loc[t,ptime])
			cntr += 1
			#~ plt.show()
		ax[0].set_title(plot_id[0]+"_"+correlation_method+"=%.2f" % corr.loc[t,ptime+"_exp_corr"])
		ax[1].set_title(plot_id[1]+"_"+correlation_method+"=%.2f" % corr.loc[t,ptime+"_ctrl_corr"])
		ax[2].set_title(plot_id[2]+"_w_Ctrl_pseudotime_"+correlation_method+"=%.2f" % ctrl_corr.loc[t,ctrl_ptime])
		
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

#~ default_output_dir = "genes_corr_w_ptime"
if not os.path.exists(output_dir):
	os.makedirs(output_dir)
	
IPython.embed()
		
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
		ctrl_ptime = raw_input("Which ctrl pseudotime would you like to correlate with? ("+ctrl_user_ptimes+ ") ")
		DEGS = pd.read_csv(DEG_path, index_col=0)
		corr["order"] = corr[ptime+"_exp_corr"].abs()
		DEGS = corr[corr.index.isin(DEGS.index)].index
		out_filename = output_dir+correlation_method+"_"+ptime+"_DEGS.pdf"
		plot_genes_of_interest(DEGS, out_filename, expression_table, annotation, pt[ptime], pt[ctrl_ptime])
		
	elif(action == "T"):
		top_n = int(raw_input("How many genes would you like to plot? "))
		ptime = raw_input("Which pseudotime would you like correlate with? ("+user_ptimes+ ") ")
		ctrl_ptime = raw_input("Which ctrl pseudotime would you like to correlate with? ("+ctrl_user_ptimes+ ") ")
		corr["order"] = corr[ptime+"_exp_corr"].abs()
		genes_of_interest = corr.sort_values(by="order", ascending=False).index[:top_n]
		out_filename = output_dir+correlation_method+"_"+ptime+".pdf"
		plot_genes_of_interest(genes_of_interest, out_filename, expression_table, annotation, pt[ptime], cpt[ctrl_ptime])
	elif(action == "I"):
		IPython.embed()

#~ if __name__ == "__main__":
	#~ main()
