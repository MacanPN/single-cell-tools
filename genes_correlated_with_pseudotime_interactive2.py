#!/usr/bin/env python 

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
from sc_pseudotime2 import *
from matplotlib.backends.backend_pdf import PdfPages
import os
import ipdb 
import pickle

from inspect import currentframe, getframeinfo

import argparse

parser = argparse.ArgumentParser(description="runs genes_correlated_with_pseudotime")
parser.add_argument("-e", "--expression-matrix", dest="expr_mat", default="/home/skevin/single_cell_tools/FACS_0407_2017_SHL_input_files/sunhye_census_matrix_20170407-comma_delimited.csv", help="gene by cell matrix of expression values", metavar="EXPR")
parser.add_argument("-c", "--cell-sets", dest="cell_sets", default="/home/skevin/single_cell_tools/FACS_0407_2017_SHL_input_files/cell_sets_0407_SHL_20180523.csv", help="cell sets", metavar="CELL_SETS")
parser.add_argument("-p", "--plot-settings", dest="plot_settings", default="/home/skevin/single_cell_tools/FACS_0407_2017_SHL_input_files/plot_setting_0407_SHL_20180212.csv", help="plot settings", metavar="PLOT_SETTINGS")
parser.add_argument("-r", "--corr-method", dest="corr_method", default="spearman", help="method of correlation (spearman or pearson)", metavar="CORR_METHOD")
parser.add_argument("-f", "--feature", dest="feature", help="feature of interest; either 'gene' or 'transcript' depending on desired output", metavar="FEATURE", required=True)
parser.add_argument("-o", "--outfile", dest="outfile", default="gene_corr_with_ptime", help="a name to give to the output file", metavar="OUTFILE")
parser.add_argument("-pt", "--pseudotime", dest="pseudotime", help="experimental cells. a list of pseudotime values. Can accept multiple values", metavar="PTIME", nargs='+', required=True)
parser.add_argument("-cpt", "--control-pseudotime", dest="ctrl_pseudotime", help="control cells. a list of pseudotime values. Can accept multiple values", metavar="PTIME", nargs='+')


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
if options.ctrl_pseudotime:
	ctrl_pseudotime_files = sorted(options.ctrl_pseudotime) 
sett = settings(settings_file, cellset_file)
expression_table, annotation = read_expression(expression_file, sett, min_expression = 0.1, min_cells = 5)
#~ expression_table = trx_expression_table.copy()

def save_obj(output_dir, obj, name):
    with open(output_dir+'/'+ name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def symbols_from_geneids(geneids):
	mg = mygene.MyGeneInfo()
	# ~ for i in enumerate(genes): 
	gene_info = mg.querymany(geneids.index, scopes='ensembl.gene', fields='symbol')
	gene_info[:] = [d for d in gene_info if d.get('notfound') != True]
	symbols = [d['symbol'] for d in gene_info]
	return(symbols)

def load_obj(output_dir, name ):
    with open(output_dir+ name + '.pkl', 'rb') as f:
        return pickle.load(f)

def list_gene_trx(output_dir, expression_table):
  gene_expression_file = output_dir+"gene_expression.csv"
  if not os.path.isfile(gene_expression_file):
  	print("creating gene expression table; saving as "+gene_expression_file)
  	gene_trx_dic = get_gene_transcript_dic(expression_table)
  	gene_expression_table =  trx_to_gene_exp_table(expression_table, gene_trx_dic)
  	gene_expression_table.to_csv(gene_expression_file, sep = "\t")
  	save_obj(output_dir, gene_trx_dic, 'gene_trx_dic')
  else:
  	gene_expression_table = pd.read_csv(gene_expression_file, index_col=0, sep ="\t")
  	gene_trx_dic = load_obj(output_dir, 'gene_trx_dic')
  return(gene_expression_table, gene_trx_dic)
  
gene_expression_table, gene_trx_dic = list_gene_trx(output_dir, expression_table)
	
def read_pts(pseudotime_files):
  pt = [read_pseudotime_from_file(i) for i in pseudotime_files]
  ptime_titles = [i.replace(".csv", "").rsplit("/")[-1] for i in pseudotime_files]
  pt = dict(zip(ptime_titles, pt))
  return pt
  
def read_corr(pt, expression_table, annotation, gene_trx_dic, cell_set_flag=None, feature=None, correlation_method=None):
  corr = [get_correlation_with_pseudotime(x, expression_table, annotation, gene_trx_dic, cell_set_flag, feature = feature, method=correlation_method) for x in pt.values()]
  corr = pd.concat(corr, axis=1)
  # 
  return corr

def load_ctrl_pseudotimes(pseudotime_files, ctrl_pseudotime_files=None, expression_table=None, annotation=None, gene_trx_dic=None, feature="transcript", correlation_method = "spearman"):
  # ipdb.set_trace()
  #load ctrl pseudotime objects if control files supplied
  if ctrl_pseudotime_files is None:
  	print("no ctrl pseudotime files supplied!")
  	ctrl_user_ptimes = "none"
  	cpt = "none"
  else:
  	# read in control pseudotime files
  	cpt = [read_pseudotime_from_file(i) for i in ctrl_pseudotime_files]
  	cptime_titles = [i.replace(".csv", "").rsplit("/")[-1] for i in ctrl_pseudotime_files]
  	cpt = dict(zip(cptime_titles, cpt))
  	if (options.feature == "transcript" or options.feature == "transcripts"):
  		ctrl_correlation_file = "_".join(cptime_titles)+"_"+correlation_method+"_transcript_correlation.csv"
  	elif (options.feature == "gene" or options.feature == "genes"):
  		ctrl_correlation_file = "_".join(cptime_titles)+"_"+correlation_method+"_symbol_correlation.csv"
  
  	#~ gene_exp_file = "_".join(cptime_titles)+"_"+correlation_method+"_gene_expression.csv"
  
  	# read correlation files from similarly named files
  	if os.path.exists(ctrl_correlation_file):
  		print(ctrl_correlation_file)
  		ctrl_corr = pd.read_csv(ctrl_correlation_file, sep="\t", index_col=0)
  
  	# check if control correlation files have already been read in
  	try:
  		ctrl_corr
  	except:
  		pt = read_pts(ctrl_pseudotime_files)
  		ctrl_corr = read_corr(pt, expression_table, annotation, gene_trx_dic, cell_set_flag="ctrl", feature=feature, correlation_method="spearman")
  		ctrl_corr.columns = sorted(cpt.keys())
  		#~ exp.to_csv(gene_exp_file, sep="\t")
  		ctrl_corr.to_csv(ctrl_correlation_file, sep="\t")
  	else:
  		if sorted(cpt.keys()) == list(ctrl_corr.columns):
  			pass
  		else:
  			print("column names do not match!")
  
  	ctrl_user_ptimes = ' '.join(cpt.keys())
  
  return(corr, cpt)

def load_exp_pseudotimes(pseudotime_files, cpt="none", expression_table=None, annotation=None, gene_trx_dic=None, feature="transcript", correlation_method = "spearman"):
  # ipdb.set_trace()
  # read in pseudotime files
  pt = [read_pseudotime_from_file(i) for i in pseudotime_files]
  ptime_titles = [i.replace(".csv", "").rsplit("/")[-1] for i in pseudotime_files]
  pt = dict(zip(ptime_titles, pt))
  
  if feature == "transcript":
  	correlation_file = "_".join(ptime_titles)+"_"+correlation_method+"_transcript_correlation.csv"
  elif feature == "gene":
  	correlation_file = "_".join(ptime_titles)+"_"+correlation_method+"_symbol_correlation.csv"

  print(correlation_file)
  
  gene_exp_file = "_".join(ptime_titles)+"_"+correlation_method+"_gene_expression.csv"
  if os.path.isfile(gene_exp_file):
  	exp = pd.read_csv(gene_exp_file, sep = "\t", index_col=0)
  
  # read correlation files from similarly named files
  if os.path.exists(correlation_file):
  	corr = pd.read_csv(correlation_file, sep="\t", index_col=0)
  
  if cpt == "none":

    pt = read_pts(pseudotime_files)
    corr = read_corr(pt, expression_table, annotation, gene_trx_dic, cell_set_flag="exp", feature="transcript", correlation_method="spearman")
    corr_columns = []
    for i in sorted(pt.keys()):
      corr_columns += [i+"_exp_corr"]
    corr.columns = corr_columns
    corr.to_csv(correlation_file, sep="\t")
  else:
    pt = read_pts(pseudotime_files)
    corr = read_corr(pt, expression_table, annotation, gene_trx_dic, cell_set_flag="mix", feature="transcript", correlation_method="spearman")
    corr_columns = []
    for i in sorted(pt.keys()):
      corr_columns += [i+"_exp_corr"]
      corr_columns += [i+"_ctrl_corr"]
      corr.columns = corr_columns
    corr.to_csv(correlation_file, sep="\t")

  corr_columns = []
  for i in sorted(pt.keys()):
  	corr_columns += [i+"_exp_corr"]
  	corr_columns += [i+"_ctrl_corr"]
  
  return(corr, pt)

corr, pt = load_exp_pseudotimes(pseudotime_files, expression_table = expression_table, feature="transcript")

## function finds genes within threshold
def genes_within_threshold(rbkd_thresh, corr):
	
	if ctrl_ptime: 
		rbkd_thresh = rbkd_thresh+"_exp_corr"
		opposite_sign = corr[(corr[rbkd_thresh]*ctrl_corr[ctrl_ptime] < 0)].index
		small_abs     = ctrl_corr[ctrl_corr[ctrl_ptime].abs() < lt].index
		good_corr_in_knockdown = corr[(corr[rbkd_thresh].abs() >= ht)].index
		genes_of_interest = corr.loc[opposite_sign.union(small_abs).intersection(good_corr_in_knockdown)]
	else: 
		rbkd_thresh = rbkd_thresh+"_exp_corr"
		good_corr_in_knockdown = corr[(corr[rbkd_thresh].abs() >= ht)].index
		genes_of_interest = corr.loc[good_corr_in_knockdown]

	return(genes_of_interest)

## function plots genes of interest (pd.Index) into pdf
def plot_genes_of_interest(genes_of_interest, out_filename, expression_table, annotation, ptime, pt, ctrl_pseudotime=None, squeeze=True):

	if ctrl_pseudotime is None:
		plot_id = ["exp"]*len(pt.keys())
	else:
		plot_id = ["exp", "Ctrl_wo_RBKD", "Ctrl_alone"]
	out_genes = out_filename.replace(".pdf", ".csv")
	og = open(out_genes, 'w')
	og.write("gene"+"\t"+"\t".join(pt.keys())+"\n")
	pp = PdfPages(out_filename)
	for i,t in enumerate(genes_of_interest):
		fig, ax = plt.subplots(1,len(plot_id), figsize=(15,5), sharey="row", squeeze=squeeze) #define common y axis for set of plots (treatments)
		print(i,t)
		title = t
		try:
			mg = mygene.MyGeneInfo()
			gene_info = mg.querymany(t, scopes='ensembl.transcript')[0]
			title = t + "  "+ gene_info["symbol"]
			title += "  ("+gene_info["name"]+")"
			correlations = [corr.loc[t,n+"_exp_corr"] for n in pt.keys()]
			correlations = "\t".join(map(str, correlations))
			
			# ~ correlations = [corr.loc[t,pt.keys()[i]] for i in len(pt.keys())]
			
			og.write(title+"\t"+correlations+"\n")
		except:
			pass
		fig.suptitle(title)
		cntr = 0
		
		if ctrl_pseudotime is None:
			while cntr < len(plot_id):
				plot_gene_with_pseudotime(expression_table, pt[list(pt)[cntr]], t, annotation, ax=ax[0][0+cntr], plot_id=plot_id[cntr])
				cntr += 1
			# ~ for i in pt.names

			for key in pt:
				ax[0][list(pt).index(key)].set_title(key+"_"+correlation_method+"=%.2f" % corr.loc[t,key+"_exp_corr"])
				# ax[0][key].set_title(key+"_"+correlation_method+"=%.2f" % corr.loc[t,key+"_exp_corr"])



		else: 
			while cntr < len(plot_id):
				plot_gene_with_pseudotime(expression_table, pt[ptime], t, annotation, ax=ax[0+cntr], plot_id=plot_id[cntr], ctrl_pseudotime=ctrl_pseudotime)
				cntr += 1
				#~ plt.show()
			ax[0].set_title(plot_id[0]+"_"+correlation_method+"=%.2f" % corr.loc[t,ptime+"_exp_corr"])
			ax[1].set_title(plot_id[1]+"_"+correlation_method+"=%.2f" % corr.loc[t,ptime+"_ctrl_corr"])
			ax[2].set_title(plot_id[2]+"_w_Ctrl_pseudotime_"+correlation_method+"=%.2f" % ctrl_corr.loc[t,ctrl_ptime])

		
		plt.tight_layout()
		plt.xlabel("pseudotime")
		plt.ylabel("log 2 expression")
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

if options.feature == "gene":
	expression_table = gene_expression_table

# main loop / choosing action
while True:
	question = """Choose from following:
	[C]	Create Correlation file from New Pseudotime Files
	[D]	Plot User-Supplied Genes
	[T]	Plot Top N Features (Genes or Transcripts) with highest correlation
	[X]	Exit
	"""
	action = input(question).upper()
	if(action == "X"):
		break
		#~ exit()
		#~ FACS_0407_2017_SHL_input_files/DEGS_day_12.csv
	elif(action == "C"):
		ptime_paths = input("provide path(s) to new pseudotime files ").split(",")
		corr_out = input("provide filename for new correlation files ")
		## block of code to calculate correlations
		pt = map(read_pseudotime_from_file, ptime_paths)
		corr = [get_correlation_with_pseudotime(x, expression_table, gene_trx_dic, method=correlation_method) for x in pt]
		corr.to_csv(corr_out+"_"+correlation_method+"_correlation.csv", sep="\t")
		ptime_titles = [i.replace(".csv", "").rsplit("/")[-1] for i in ptime_paths]
		ptimes = dict(zip(ptime_titles, ptime_paths))
		user_ptimes = ' '.join(ptime_titles)
		corr.columns = ptime_titles
		
	elif(action == "D"):
		DEG_path = input("provide path to differentially expressed genes ")
		ptime = input("Which pseudotime would you like correlate with? ("+user_ptimes+ ") ")
		ctrl_ptime = input("Which ctrl pseudotime would you like to correlate with? ("+ctrl_user_ptimes+ ") ")
		DEGS = pd.read_csv(DEG_path, index_col=0, header=None)
		if options.feature == "gene":
			DEGS = symbols_from_geneids(DEGS)
		elif options.feature == "transcript":
			DEGS = DEGS.index
		corr["order"] = corr[ptime+"_exp_corr"].abs()

		DEGS = corr[corr.index.isin(DEGS)].index
		out_filename = output_dir+correlation_method+"_"+ptime+"_DEGS.pdf"

		if ctrl_ptime == '':
			if len(pt) == 1:

			  plot_genes_of_interest(DEGS, out_filename, expression_table, annotation, ptime, pt, squeeze=False)
        # plot transcripts of interest
				#~ plot_genes_of_interest(DEGS, out_filename, expression_table, annotation, ptime, pt, squeeze=False)
				
			else:
				plot_genes_of_interest(DEGS, out_filename, expression_table, annotation, ptime, pt)
				# plot transcripts of interest
				#~ plot_genes_of_interest(DEGS, out_filename, expression_table, annotation, ptime, pt)
		else:
			plot_genes_of_interest(DEGS, out_filename, expression_table, annotation, ptime, pt, cpt[ctrl_ptime])
			# plot transcripts of interest
			#~ plot_genes_of_interest(DEGS, out_filename, expression_table, annotation, ptime, pt, cpt[ctrl_ptime])
		
	elif(action == "T"):
		top_n = int(input("How many genes would you like to plot? "))
		ptime = input("Which pseudotime would you like to order by? ("+user_ptimes+ ") ")
		ctrl_ptime = input("Which ctrl pseudotime would you like to correlate with? Leave blank if no control present. ("+ctrl_user_ptimes+ ") ")

		ht = input("set upper threshold (def. 0.3) ")
		if not ht == '':
		  ht = float(ht)
		else:
		  ht = 0.3
		lt = float(input("set lower threshold (def. 0.2) Leave blank if no control present. ") or 0)
		threshold_set = input("retain genes above upper threshold in ("+user_ptimes+ ") (split with , for mult.) ").split(",")
		corr["order"] = (corr[ptime+"_exp_corr"]).abs()
		# ~ corr["order"] = (corr[ptime]).abs()
	
		if len(threshold_set) > 1:

			corrs = map(genes_within_threshold, threshold_set, corr)
			genes_of_interest = corr.loc[pd.concat(corrs, axis=1, join='inner').index]

		else:
			genes_of_interest = genes_within_threshold(threshold_set[0], corr)

		if top_n > len(genes_of_interest):
			print("error! number of genes requested exceeds number of genes matching filtering criteria ("+str(len(genes_of_interest))+")")
			pass
		elif ctrl_ptime == '':
			genes_of_interest = genes_of_interest.sort_values(by="order", ascending=False).index[:top_n]
			out_filename = output_dir+correlation_method+"_"+ptime+"_top_"+str(top_n)+"_genes.pdf"
			if len(pt) == 1:
				plot_genes_of_interest(genes_of_interest, out_filename, expression_table, annotation, ptime, pt, squeeze=False)
				
				# plot transcripts of interest
				#~ plot_genes_of_interest(genes_of_interest, out_filename, expression_table, annotation, ptime, pt, squeeze=False)

			else:
				plot_genes_of_interest(genes_of_interest, out_filename, expression_table, annotation, ptime, pt)
				# plot transcripts of interest
				#~ plot_genes_of_interest(genes_of_interest, out_filename, expression_table, annotation, ptime, pt)
		else:
			genes_of_interest = genes_of_interest.sort_values(by="order", ascending=False).index[:top_n]
			out_filename = output_dir+correlation_method+"_"+ptime+"_top_"+str(top_n)+"_genes.pdf"
			plot_genes_of_interest(genes_of_interest, out_filename, expression_table, annotation, ptime, pt, cpt[ctrl_ptime])
			# plot transcripts of interest
			#~ plot_genes_of_interest(genes_of_interest, out_filename, expression_table, annotation, ptime, pt, cpt[ctrl_ptime])
	elif(action == "I"):
		IPython.embed()

#~ if __name__ == "__main__":
	#~ main()