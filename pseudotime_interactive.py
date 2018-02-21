#!/usr/bin/python

import pandas as pd
import matplotlib.pyplot as plt
import sys
import numpy as np
import scipy as sc
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
	
	scipy_linkage_methods = [ "complete", "average", "single", "centroid", "median", "ward"]
	question = "Which PCs would you like to use for clustering? [type comma separated list, list can also include ranges 1-5,8] "
	cluster_on_pcs = list_from_ranges(raw_input(question))
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
		plot_hierarchical_clustering(subset_PC_expression[cluster_on_pcs], subset_annotation, method=method)
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
		plot_hierarchical_clustering(PC_expression[cluster_on_pcs], annotation, method=method)
	
	return clusters

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
	
	

# main loop / choosing action
while True:
	question = """Choose from following:
	[H]	Plot Hierarchical Clustering
	[P]	Plot PCA
	[L]	Assign time clusters according to time Labels (like day_4 ... )
	[C]	Assign time clusters using hierarchical clustering
	[D]	Find Most Correlated PCs
	[N]	Normalize centroids
	[S]	Calculate pseudotime for cells using times assigned to clusters
	[O]	Output clusters (so they can be copied to a file)
	[F]	Save generated pseudotime to file
	[X]	Exit
	"""
	action = raw_input(question).upper()
	if(action == "X"):
		break
		#~ exit()
	elif(action == "H"):
		print("plotting...\n dendrogram will be saved as a .pdf shortly")
		plot_all_hierarchical_clusterings(PC_expression, annotation, sett)
	elif(action == "P"):
		colnm, colvalp = retrieve_subset_param()
		pcs = map(int,raw_input("Which PCs would you like on the plot? (type comma separated list, such as 1,3,4) ").split(","))
		sett.pcs = pcs
		print("plotting...\n the plot will open in your web browser shortly")
		if not all(colvalp):
			plot_3d_pca(PC_expression, annotation, sett, clusters = clusters)
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
		pcs = map(int,raw_input("Which PCs would you like on the plot? (type comma separated list, such as 1,3,4) ").split(","))
		find_pseudotime(subset_PC_expression, subset_annotation, pca, sett, pcs)
		print("Time clusters were assigned according to labels")
		
	elif(action == "C"):
		colnm, colval = retrieve_subset_param()
		subset_annotation, subset_PC_expression = subset_pc_expression(PC_expression, colnm, colval)
		subset_clusters = assign_time_clusters_using_clustering(colnm, colval)
		print("Time clusters were assigned according to hierarchical clustering")
		plt.savefig("hierarchical_clustering.pdf")
		
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
		
	elif(action=="F"):
		if("pseudotime" not in globals()):
			print("Pseudotime was not yet generated!")
			continue
		filename = raw_input("Enter file name you'd like to save pseudotime as (preferably ending with .csv) ")
		pseudotime.to_csv(filename, sep="\t")
	elif(action == "I"):
		IPython.embed()
	



#HTML_pal = ['#a50026','#d73027','#f46d43','#fdae61','#fee08b','#ffffbf','#d9ef8b','#a6d96a','#66bd63','#1a9850','#006837']
#IPython.embed()
#annotation["color"] = [HTML_pal[i] for i in color_indices]
