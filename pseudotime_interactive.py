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

from sc_pseudotime import *

expression_file = sys.argv[1]
cellset_file    = sys.argv[2]
settings_file   = sys.argv[3]
n_pca = 20
# read settings and cell_set files
sett = settings(settings_file, cellset_file)
# read expression table
expression_table, annotation = read_expression(expression_file, sett)
# calculate PCA
PC_expression,pca = run_PCA(expression_table, annotation, n_pca)

clusters = None
annotation["name"] = "day "+annotation["day"].astype(str)

def assign_time_clusters_using_clustering():
	scipy_linkage_methods = [ "complete", "average", "single", "centroid", "median", "ward"]
	question = "Which PCs would you like to use for clustering? [type comma separated list, list can also include ranges 1-5,8] "
	cluster_on_pcs = list_from_ranges(raw_input(question))
	number_of_clusters = int(raw_input("How many clusters would you like to generate? "))
	method = raw_input("Which clustering would you like to use: "+", ".join(scipy_linkage_methods)+": ")
	if method not in scipy_linkage_methods:
		print("clustering method not supported (spelling error?)")
		return
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
	return clusters

def print_clusters(clusters):
	if(clusters == None):
		print("Clusters were not identified yet!")
		return
	print("Follows list of clusters: (first column is the time assigned to the cluster)")
	for cl in clusters:
		print(str(cl[0])+"\t"+"\t".join(cl[1]))

# main loop / choosing action
while True:
	question = """Choose from following:
	[P]	Plot PCA
	[L]	Assign time clusters according to time Labels (like day_4 ... )
	[C]	Assign time clusters using hierarchical clustering
	[S]	Calculate pseudotime for cells using times assigned to clusters
	[O]	Output clusters (so they can be copied to a file)
	[F]	Save generated pseudotime to file
	[X]	Exit
	"""
	action = raw_input(question).upper()
	if(action == "X"):
		break
	elif(action == "P"):
		pcs = map(int,raw_input("Which PCs would you like on the plot? (type comma separated list, such as 1,3,4) ").split(","))
		sett.pcs = pcs
		print("plotting...\n the plot will open in your web browser shortly")
		plot_3d_pca(PC_expression, annotation, sett, clusters = clusters)
	elif(action == "L"):
		clusters = time_clusters_from_annotations(annotation)
		print("Time clusters were assigned according to labels")
	elif(action == "C"):
		clusters = assign_time_clusters_using_clustering()
		print("Time clusters were assigned according to hierarchycal clustering")
	elif(action == "S"):
		pseudotime = calcutale_pseudotime_using_cluster_times(PC_expression, annotation, clusters, sett)
	elif(action == "O"):
		print_clusters(clusters)
	elif(action=="F"):
		if("pseudotime" not in globals()):
			print("Pseudotime was not yet generated!")
			continue
		filename = raw_input("Enter file name you'd like to save pseudotime as (preferably ending with .csv)")
		pseudotime.to_csv(filename, sep="\t")
	elif(action == "I"):
		IPython.embed()
	


#HTML_pal = ['#a50026','#d73027','#f46d43','#fdae61','#fee08b','#ffffbf','#d9ef8b','#a6d96a','#66bd63','#1a9850','#006837']
#IPython.embed()
#annotation["color"] = [HTML_pal[i] for i in color_indices]
