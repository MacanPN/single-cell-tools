#!/usr/bin/python

import pandas as pd
import matplotlib
matplotlib.use('Agg') #for handling display variable on mac (SHL!)
import matplotlib.pyplot as plt
import sys
import numpy as np
import scipy as sc
import scipy.cluster.hierarchy as sch
import sklearn
import seaborn as sns
import math
import IPython
import argparse
import os
import plotly
#~ import ipdb

from sc_pseudotime import *

parser = argparse.ArgumentParser(description="runs pseudotime_interactive")
parser.add_argument("-e", "--expression-matrix", dest="expr_mat", default="~/single_cell_tools/example_input_files/transcripts.tpm_census_matrix-comma-delimited.csv", help="gene by cell matrix of expression values", metavar="EXPR")
parser.add_argument("-c", "--cell-sets", dest="cell_sets", default="~/single_cell_tools/example_input_files/cell_sets.csv", help="cell sets", metavar="CELL_SETS")
parser.add_argument("-p", "--plot-settings", dest="plot_settings", default="~/single_cell_tools/example_input_files/plot_settings.csv", help="plot settings", metavar="PLOT_SETTINGS")
parser.add_argument("-n", "--session-name", dest="session_name", help="a name to give to this analysis session for reproducbility", metavar="SESSION_NAME", required=True)


try:
    options = parser.parse_args()
except SystemExit as err: 
	if err.code == 2: 
		parser.print_help()
		sys.exit(0)
 
expression_file = os.path.expanduser(options.expr_mat)
cellset_file    = os.path.expanduser(options.cell_sets)
settings_file   = os.path.expanduser(options.plot_settings)
n_pca = 20
# read settings and cell_set files
sett = settings(settings_file, cellset_file)
# read expression table
expression_table, annotation = read_expression(expression_file, sett)
# calculate PCA

PC_expression,pca = run_PCA(expression_table, annotation, n_pca)

clusters = None
annotation["name"] = "day "+annotation["day"].astype(str)

def assign_time_clusters_using_clustering(colnm=None, colval=None):
	pc_set = "Which PCs would you like to use for clustering? [type comma separated list, list can also include ranges 1-5,8] "
	scipy_linkage_methods = [ "complete", "average", "single", "centroid", "median", "ward"]
	cluster_on_pcs = list_from_ranges(raw_input(pc_set))
	number_of_clusters = int(raw_input("How many clusters would you like to generate? "))
	method = raw_input("Which clustering would you like to use: "+", ".join(scipy_linkage_methods)+": ")
	if method not in scipy_linkage_methods:
		print("clustering method not supported (spelling error?)")
		return
	if '' not in colval:
		linkage = sc.cluster.hierarchy.linkage(subset_PC_expression[cluster_on_pcs], method=method)
		clusters_without_time = get_cluster_labels(linkage, number_of_clusters, subset_PC_expression.index)
		cluster_colors = ["blue", "red", "orange", "purple", "green", "brown", "black", "gray", "lawngreen", "magenta", "lightpink", "indigo", "lightblue", "lightgoldenrod1", "mediumpurple2"]
		print("Now plotting clusters")
		change_annotation_colors_to_clusters(clusters_without_time, subset_annotation, cluster_colors)
		clusters = []
		sett.pcs = cluster_on_pcs[:3]
		plot_3d_pca(subset_PC_expression, subset_annotation, sett)
		for i in range(0,number_of_clusters):
			time = float(raw_input("Assign time for cluster shown in "+cluster_colors[i]+": "))
			clusters.append( (time,subset_annotation.loc[subset_annotation["color"]==cluster_colors[i]].index) )
		clusters.sort(key=lambda by_first: by_first[0])
		dendro = plot_hierarchical_clustering(subset_PC_expression[cluster_on_pcs], subset_annotation, method=method)
	else:
		linkage = sc.cluster.hierarchy.linkage(PC_expression[cluster_on_pcs], method=method)
		clusters_without_time = get_cluster_labels(linkage, number_of_clusters, PC_expression.index)
		cluster_colors = ["blue", "red", "orange", "purple", "green", "brown", "black", "gray", "lawngreen", "magenta", "lightpink", "indigo", "lightblue", "lightgoldenrod1", "mediumpurple2"]
		print("Now plotting clusters")
		change_annotation_colors_to_clusters(clusters_without_time, annotation, cluster_colors)
		clusters = []
		sett.pcs = cluster_on_pcs[:3]
		plot_3d_pca(PC_expression, annotation, sett)
		for i in range(0,number_of_clusters):
			time = float(raw_input("Assign time for cluster shown in "+cluster_colors[i]+": "))
			clusters.append( (time,annotation.loc[annotation["color"]==cluster_colors[i]].index) )
		clusters.sort(key=lambda by_first: by_first[0])
		dendro = plot_hierarchical_clustering(PC_expression[cluster_on_pcs], annotation, method=method)
	
	return clusters, dendro

def print_clusters(clusters):
	if(clusters == None):
		print("Clusters were not identified yet!")
		return
	print("Follows list of clusters: (first column is the time assigned to the cluster)")
	for cl in clusters:
		print(str(cl[0])+"\t"+"\t".join(cl[1]))
	# alternative cluster_print method
	#~ clusters_df = pd.DataFrame()
	#~ for i, (a, b) in enumerate(clusters):
		#~ enum_df = pd.DataFrame(clusters[i][1])
		#~ enum_df.columns = [str(a)]
		#~ clusters_df = pd.concat([clusters_df, enum_df], axis=1)

def retrieve_subset_param():
	colnm = raw_input("What metadata should be used to subset the data? (ex. treatment, age, etc.) ")
	colval = raw_input("What values should be used to subset the data? (ex. shCtrl, sh842,). Providing no value will prevent subsetting ").split(",")
	return colnm, colval

def subset_pc_expression(pc_expression, colnm, colval):
	if not all(colval):
		# ~ IPython.embed()
		# ~ clusters_without_time = get_cluster_labels(linkage, number_of_clusters, subset_PC_expression.index)
		# ~ IPython.embed()
		# ~ cluster_colors = ["blue", "red", "orange", "purple", "green", "brown", "black", "gray", "lawngreen", "magenta", "lightpink", "indigo", "lightblue", "lightgoldenrod1", "mediumpurple2"]
		# ~ change_annotation_colors_to_clusters(clusters_without_time, subset_annotation, cluster_colors)
		return annotation, PC_expression
	else:
		if colnm not in annotation.columns:
			print("metadat not recognized (spelling error?)")
			colnm, colval = retrieve_subset_param()
		subset_annotation = annotation[annotation[colnm].isin(colval)]
		# add day0 cells to all subset_annotations if not removed in plot_settings
		day0_annotation = annotation[annotation["day"]==0.0]
		excluded_day0 = day0_annotation[-day0_annotation.isin(subset_annotation)].dropna()
		if (not day0_annotation.empty):
			subset_annotation = subset_annotation.append(excluded_day0)
		subset_PC_expression = PC_expression.loc[subset_annotation.index.values]
		return subset_annotation, subset_PC_expression

def find_discrim_pcs(subset_pc_expression, annotation):
	
	shRBKD_cells = annotation.loc[annotation['treatment'] != "shCtrl"].index
	vals_shRBKD = PC_expression.loc[shRBKD_cells,:]
	
	shCtrl_cells = annotation.loc[annotation['treatment'] == "shCtrl"].index
	vals_shCtrl = PC_expression.loc[shCtrl_cells,:]
	
	for i in range(1, 20):
		test = (abs(vals_shCtrl.loc[:,i].mean()-vals_shRBKD.loc[:,i].mean()))/np.std(vals_shCtrl.loc[:,i].append(vals_shRBKD.loc[:,i]))
		print(test)

def normalize_centroids(subset_pc_expression):
	print("provide control group: ")
	ctrl_colnm, ctrl_colval = retrieve_subset_param()
	pcs = map(int,raw_input("Which PCs would you like on the plot? (type comma separated list, such as 1,3,4) ").split(","))
	sett.pcs = pcs
	
	ctrl_annotation, ctrl_pc_expression = subset_pc_expression(PC_expression, colnm, colval)
	ctrl_clusters = time_clusters_from_annotations(ctrl_annotation)
	ctrl_cntrds = get_cluster_centroids(ctrl_pc_expression, ctrl_clusters)
	try:
		subset_PC_expression
	except:
		print("error")
	else:
		fig = plot_3d_pca(subset_PC_expression, subset_annotation, sett, clusters = subset_clusters)
		trace = fig["data"][-1]
		sub_ctrl_cntrds = ctrl_cntrds.iloc[sett.pcs,:]
		cntrds = pd.DataFrame([])
		coords = ["x", "y", "z"]
		for i,j in zip(coords, sub_ctrl_cntrds.index):
			trace[i] = trace[i] - (sub_ctrl_cntrds.loc[j,] - sub_ctrl_cntrds.loc[j,0])
		fig["data"][-1] = trace
	#~ IPython.embed()
	return(fig)
	#~ trace = record_centroids(clusters, comb, settings)
	#~ centroids = pd.DataFrame([])
	#~ coords = ["x", "y", "z"]
	#~ for i in coords:
		#~ centroids = pd.append([trace[i])

## save genes correlated with a PC to file
def get_isoforms_correlated_with_pc(pc, n):
	pc = int(pc)
	df = pd.Series(pca.components_[pc], index=expression_table.columns)
	# ~ df = pd.DataFrame(pca.components_[pc], pc, index=expression_table.columns)
	df = df.reindex(df.abs().sort_values(inplace=False, ascending=False).index).iloc[0:n]
	return df
	
def get_isoforms_correlated_pc_set(pca, expression_table, pcs, n, filename):
	data = [get_isoforms_correlated_with_pc(i, n) for i in pcs]
	data = pd.concat(data)
	csv_filename = filename+"_pcs_"+"_".join(map(str, pcs))+".csv"
	data.to_csv(csv_filename, sep=",")
	return data
	

# main loop / choosing action
while True:
	question = """Choose from following:
	[H]	Plot Hierarchical Clustering
	[P]	Plot PCA
	[L]	Assign time clusters according to time Labels (like day_4 ... )
	[C]	Assign time clusters using hierarchical clustering
	[D]	Find Most Correlated and Most Discriminating (treat v ctrl) PCs
	[N]	Normalize centroids
	[G]	Plot 2d PCA of Marker Genes Colored by Expression
	[S]	Calculate pseudotime for cells using times assigned to clusters
	[O]	Output clusters (so they can be copied to a file)
	[F]	Save generated pseudotime to file
	[M]	Save features correlated with pc to file
	[X]	Exit
	"""
	action = raw_input(question).upper()
	if(action == "X"):
		break
		#~ exit()
	elif(action == "H"):
		number_of_clusters = int(raw_input("How many clusters would you like to generate? "))
		sett.num_clusters = number_of_clusters
		print("plotting...\n dendrogram will be saved as a .pdf shortly")
		try:
			subset_annotation
		except:
			plot_all_hierarchical_clusterings(PC_expression, annotation, sett)
		else:
			plot_all_hierarchical_clusterings(PC_expression, subset_annotation, sett)
	elif(action == "P"):
		colnm, colvalp = retrieve_subset_param()
		pcs = map(int,raw_input("Which PCs would you like on the plot? (type comma separated list, such as 1,3,4) ").split(","))
		sett.pcs = pcs
		print("plotting...\n the plot will open in your web browser shortly")
		if not all(colvalp):
			fig = plot_3d_pca(PC_expression, annotation, sett, clusters = clusters)
		else:
			try:
				subset_PC_expression
			except:
				subset_annotation, subset_PC_expression = subset_pc_expression(PC_expression, colnm, colvalp)
				plot_3d_pca(subset_PC_expression, subset_annotation, sett, clusters = clusters)
				del subset_annotation, subset_PC_expression
			else:
				if (colvalp == colval):
					plot_3d_pca(subset_PC_expression, subset_annotation, sett, clusters = subset_clusters)
				elif set(colvalp).issubset(colval):
					#~ ipdb.set_trace()
					old_colors = subset_annotation["color"]
					subset_annotation, subset_PC_expression = subset_pc_expression(subset_PC_expression, colnm, colvalp)
					subset_annotation.loc[:,"color"] = old_colors[old_colors.index.isin(subset_annotation.index)]
					#~ IPython.embed()
					new_clusters = [(i, c[c.isin(subset_annotation.index)]) for i,c in subset_clusters]
					plot_3d_pca(subset_PC_expression, subset_annotation, sett, clusters = new_clusters)
					del subset_annotation, subset_PC_expression
				else:
					subset_annotation, subset_PC_expression = subset_pc_expression(PC_expression, colnm, colvalp)
					plot_3d_pca(subset_PC_expression, subset_annotation, sett, clusters = clusters)
					del subset_annotation, subset_PC_expression
	elif(action == "L"):
		colnm, colval = retrieve_subset_param()
		subset_annotation, subset_PC_expression = subset_pc_expression(PC_expression, colnm, colval)
		subset_clusters = time_clusters_from_annotations(subset_annotation)
		pcs = map(int,raw_input("Which PCs would you like on the plot? (type comma separated list, such as 1,3,4) ").split(","))
		print("Time clusters were assigned according to labels")
	
	elif(action == "D"):
		colnm, colval = retrieve_subset_param()
		subset_annotation, subset_PC_expression = subset_pc_expression(PC_expression, colnm, colval)
		# ~ pcs = map(int,raw_input("Which PCs would you like on the plot? (type comma separated list, such as 1,3,4) ").split(","))
		find_pseudotime(subset_PC_expression, subset_annotation, pca, sett)
		print("Showing PCS most correlated with time")
		
	elif(action == "C"):
		colnm, colval = retrieve_subset_param()
		subset_annotation, subset_PC_expression = subset_pc_expression(PC_expression, colnm, colval)
		subset_clusters, dendro = assign_time_clusters_using_clustering(colnm, colval)
		print("Time clusters were assigned according to hierarchical clustering")
		filename = raw_input("Enter file name you'd like to save clustering plot as (preferably ending with .pdf) ")
		plt.savefig(filename)
		
	elif(action == "N"):
		test = normalize_centroids(subset_pc_expression)
		url = plotly.offline.plot(test, filename="normalize_centroids.html", validate=False, auto_open=False)
		print(url)
		
	elif(action == "S"):
		colnm, colval = retrieve_subset_param()
		subset_annotation, subset_PC_expression = subset_pc_expression(PC_expression, colnm, colval)
		pseudotime, centroids = calculate_pseudotime_using_cluster_times(subset_PC_expression, subset_annotation, subset_clusters, sett)
		
	elif(action == "O"):
		print_clusters(subset_clusters)
		
	elif(action == "G"):
		colnm, colval = retrieve_subset_param()
		subset_annotation, subset_PC_expression = subset_pc_expression(PC_expression, colnm, colval)
		pcs = map(int,raw_input("Which PCs would you to correlate with? (type comma separated list, such as 1,3,4) ").split(","))
		sett.pcs = pcs
		bins = int(raw_input("How many bins would you like to quantile? "))
		sett.bins = bins
		marker_genes = raw_input("Which marker genes would you like to plot (type comma separated list, such as RB1,RXRG,ARR3) ").split(",")
		bin_col_dict = {}
		for i in range(0,bins):
			color = raw_input("Assign color for bin "+str(i)+": ")
			bin_col_dict.update({i:color})
		for i in marker_genes:
			plot_3d_pca(subset_PC_expression, subset_annotation, sett, clusters = clusters, genes=i, bin_col_dict=bin_col_dict, expression_table=expression_table)
		# ~ plot_marker_gene_quantile(expression_table, subset_PC_expression, subset_annotation, pcs, sett, marker_genes)
		
		
	elif(action == "M"):
		filename = raw_input("Enter file name you'd like to save correlated features to: ")
		top_n = int(raw_input("How many genes from each pc would you like to plot? "))
		pcs = map(int,raw_input("Which PCs would you to correlate with? (type comma separated list, such as 1,3,4) ").split(","))
		pc_corr_trs = get_isoforms_correlated_pc_set(pca, expression_table, pcs, top_n, filename)
		
	elif(action=="F"):
		if("pseudotime" not in globals()):
			print("Pseudotime was not yet generated!")
			continue
		filename = raw_input("Enter file name you'd like to save pseudotime as (preferably ending with .csv) ")
		pseudotime.to_csv(filename, sep="\t")
		
	elif(action=="Q"):
		cluster_dir = raw_input("Enter location of diffex csvs ")
		diffex_csvs = read_in_diffex(cluster_dir)
		plot_heatmap(expression_table, annotation, dendro)
		
	elif(action == "I"):
		IPython.embed()
	




#HTML_pal = ['#a50026','#d73027','#f46d43','#fdae61','#fee08b','#ffffbf','#d9ef8b','#a6d96a','#66bd63','#1a9850','#006837']
#IPython.embed()
#annotation["color"] = [HTML_pal[i] for i in color_indices]
