#!/usr/bin/python

## @program sc-analysis
# file sc-analysis.py
# this package is built to enable temporal analysis of single cell rna-seq data

import pandas as pd
import matplotlib.pyplot as plt
import sys
import numpy as np
import scipy as sc
import scipy.cluster.hierarchy as sch
import sklearn
from sklearn import decomposition
from sklearn import cluster
import seaborn as sns
import math
import IPython
import mygene
import statsmodels.api as sm
import copy # needed to copy class settings
from matplotlib.backends.backend_pdf import PdfPages
import plotly.plotly as py
import plotly
import plotly.graph_objs as go
#~ import plotly.graph_objs as go
#modules for heatmap plotting 
from plotly.graph_objs import *
import plotly.figure_factory as FF
from scipy.spatial.distance import pdist, squareform
import re
import matplotlib.cm as cm
import pickle

import ipdb
import os
## what modes can be script run in
run_modes = ["2d-pca-multiplot", "2d-pca-single", "3d-pca", "hierarchy", "pseudotime", "3d-pca-colored-by-clustering", "test"]
default_shape = "o"
default_day = 1.0
default_color = "gray"
time_group_prefix = "day_"
treatment_group_prefix = "sh"
cluster_group_prefix= "cluster_"
## sets with parameter look like:
# operation	set_name	parameter
# for ex.: color	day_4	clue
accepted_sets_with_parameter = ["color", "outline-color", "size", "name", "shape"]

## sets without parameter look like:
# operation	set_name
# for ex.: remove	low_read_count_cells
accepted_sets_without_parameter = ["remove", "superimpose", "superimpose-for-spearman"]

## parameters supposed to be set once
accepted_parameters = ["number_of_genes"]

def RGBToHTMLColor(rgb_tuple):
	""" convert an (R, G, B) tuple to #RRGGBB """
	hexcolor = '#%02x%02x%02x' % rgb_tuple
	# that's it! '%02x' means zero-padded, 2-digit hex values
	return hexcolor

## this class reads settings of the run and keeps them in its attributes
#  if settings file in incorrect, the init method prints error and terminates application
class settings:
	## read sets of cells later used to refer to them (to remove/color/superimpose etc...)
	#
	# format:
	#
	# cell_set_name <tab> cell <tab> cell ...
	#
	# cell_set_name <tab> cell <tab> cell ...
	def read_cell_sets(self,cellset_file):
		cell_sets = {}
		with open(cellset_file, 'rU') as f:
			for line in f:
				x = line.rstrip().split("\t")
				cell_sets[x[0]] = x[1:]
		return cell_sets
	
	def __init__(self, settings_file, cellset_file):
		self.cell_sets = self.read_cell_sets(cellset_file)
		self.sets = {}
		# inniciate all sets to empty
		for i in accepted_sets_with_parameter+accepted_sets_without_parameter:
			self.sets[i] = [] #set()
		# in some cases we'll want to keep information which PCs to plot
		self.pcs = []
		# how many genes per PC do we want to save (for gene onthology analysis)
		self.parameters = {}
		self.parameters["number_of_genes"] = "1000"
		with open(settings_file) as f:
			# first line defines name of the output
			self.result_filename = f.readline().rstrip()
			# second line defines analysis to run
			mode_line = f.readline().rstrip()
			self.run_mode = mode_line.split("\t")[0]
			# third line defines dimensions of pca plot
			dim_line = f.readline().rstrip()
			self.plot_dim = map(int, dim_line.split(","))
			if self.run_mode not in run_modes:
				print "Unkown run mode (line 2 in settings file): ",self.run_mode
				raise ValueError
			# if we're plotting pca, we want list of PCs to use
			if self.run_mode in ["2d-pca-single", "3d-pca","3d-pca-colored-by-clustering","test"]:
				self.pcs = map(int,mode_line.split("\t")[1].split(","))
				if not(
					((self.run_mode == "2d-pca-single")and(len(self.pcs)==2))
					or((self.run_mode in ["3d-pca","3d-pca-colored-by-clustering","test"])and(len(self.pcs)==3))
					):
					print "Invalid number of PCs given! ",mode_line
					raise ValueError
			# if we want to color cells by clustering, we need to grab number of clusters required
			if self.run_mode == "3d-pca-colored-by-clustering":
				self.n_clusters = int(mode_line.split("\t")[2])
				self.clustering_method = mode_line.split("\t")[3]
			# from third line onwards, the script reads different operations carried out on defined cell sets
			for line in f:
				if(line.startswith("#")):
					continue
				x = line.rstrip().split("\t")
				if(x[0] in accepted_sets_without_parameter):
					self.sets[x[0]].append(x[1]) #add(x[1])
				elif(x[0] in accepted_sets_with_parameter):
					self.sets[x[0]].append(tuple(x[1:3]))
				elif(x[0] in accepted_parameters):
					self.parameters[x[0]] = x[1]
				else:
					print "Unknown option:",line
					exit(1)
	def __copy__(self):
		return copy.deepcopy(self)


## function takes expression file and settings object and returns:
# - pd.DataFrame with [log transformed] expression values [genes expressed over min_expression in at least min_cells]
# - pd.DataFrame with annotations for each cell. Expression table and annotation table have the same rows
def read_expression(expression_file, settings, min_expression = 0.1, min_cells = 10, log_transform = True):
	# read expression
	expression_table = pd.read_csv(expression_file, sep=",", index_col = 0).transpose()
	print "Read expression table with shape:",expression_table.shape
	
	# remove genes with less then min_cells expressing them
	expressed_genes = (expression_table > min_expression).sum() > min_cells
	expression_table = expression_table.loc[ : , expressed_genes]
	print "Removed genes that are not expressed >",min_expression," in at least",min_cells ,"cells"
	print "Expression table has now shape:",expression_table.shape
	
	# remove unwanted cells
	for s in settings.sets["remove"]:
		print "Removed cells from set:",s,settings.cell_sets[s]
		expression_table.drop(settings.cell_sets[s], inplace=True, errors="ignore")

	# log transform
	if(log_transform):
		expression_table += 1
		expression_table = expression_table.apply(np.log2)
		print "Log transformed data"
	
	# create annotation table and populate it with default values
	annotation = pd.DataFrame(index=expression_table.index)
	annotation["color"] = default_color
	annotation["superimpose"] = False
	annotation["superimpose-for-spearman"] = False
	annotation["size"] = 5.0
	annotation["name"] = ""
	annotation["outline-color"] = None
	annotation["day"]=default_day
	annotation["shape"] = default_shape
	annotation["treatment"]= "none"
	annotation["cluster"]= "none"
	
	# annotate superimposed cells
	for s in settings.sets["superimpose"]:
		print "Superimposing cells from set:",s,settings.cell_sets[s]
		for i in settings.cell_sets[s]:
			annotation.loc[i,"superimpose"] = "True"
		

	
	for s in accepted_sets_with_parameter: # iterating over dictionary operation->set
		for i in settings.sets[s]: # iterating over set
			subset = set(settings.cell_sets[i[0]]).intersection(annotation.index)
			annotation.loc[subset,s] = i[1]
	
	annotation["size"] = pd.to_numeric(annotation["size"])
	# where outline color was not defined, set it to the color of the cell
	annotation.loc[annotation["outline-color"]!=annotation["outline-color"], "outline-color"] = annotation["color"]
	
	# define day, treatment, and cluster columns
	day_labels = [d for d in settings.cell_sets if d.startswith(time_group_prefix)]
	treatment_labels = [t for t in settings.cell_sets if t.startswith(treatment_group_prefix)]
	cluster_labels = [c for c in settings.cell_sets if c.startswith(cluster_group_prefix)]
	for i in day_labels:
		subset = set(settings.cell_sets[i]).intersection(annotation.index)
		annotation.loc[subset,"day"]=int(i.split("_")[1])
	for i in treatment_labels:
		subset = set(settings.cell_sets[i]).intersection(annotation.index)
		annotation.loc[subset,"treatment"]=i
	for i in cluster_labels:
		subset = set(settings.cell_sets[i]).intersection(annotation.index)
		annotation.loc[subset,"cluster"]=i.split("_")[1]
	
	# crop annotation dataframe to only rows, that are in expression table
	annotation = annotation.loc[expression_table.index]
	return (expression_table, annotation)

## runs PCA and returns:
# - PCA transformed coordinates
# - sklearn.decomposition.pca object
def run_PCA(expression_table, annotation, n_components):
	pca = decomposition.PCA(n_components=n_components, svd_solver="full")
	expression_table_for_PCA = expression_table.loc[annotation[annotation["superimpose"]==False].index]
	#~ ipdb.set_trace()
	print "Calculating PCA on table of shape:",expression_table_for_PCA.shape
	pca.fit(expression_table_for_PCA)
	print "Explained variance: ", pca.explained_variance_
	print "Explained variance ratio: ", pca.explained_variance_ratio_
	# transform expression using PCA vectors
	transformed_expression = pd.DataFrame(pca.transform(expression_table), index=expression_table.index, columns = range(1,n_components+1))
	return transformed_expression, pca

## save genes correlated with a PC to file
# ~ def get_isoforms_correlated_with_pc(pc):
	# ~ pc = int(pc)
	# ~ df = pd.Series(pca.components_[pc], index=expression_table.columns)
	# ~ df = df.reindex(df.abs().sort_values(inplace=False, ascending=False).index).iloc[0:n]
	# ~ csv_filename = filename+"_PC"+str(pc)+".csv"
	# ~ df.to_csv(csv_filename, sep="\t")
	# ~ return df
	
# ~ def get_isoforms_correlated_pc_set(pca, expression_table, pcs, n, filename):
	# ~ data = [get_isoforms_correlated_with_pc(i) for i in pcs]
	# ~ pd.concat(data)
	# ~ names = map(str, pcs)
	# ~ data = pd.DataFrame.from_items(zip(names, pca.components_[pcs]))
	# ~ data.index = expression_table.columns
	# ~ data = data.reindex(data.abs().sort_values(inplace=False, ascending=False).index).iloc[0:n]
	# ~ csv_filename = filename+"_PC"+".csv"
	# ~ data.to_csv(csv_filename, sep="\t")
	# ~ return data

## create annotation label for a point on axis if it's far enough from other points
# used internally by plotting functions
def annotate_df(row,df,min_dist,ax):
		dist = (df - row).abs().sum(axis=1).sort_values()[1]
		if(dist > min_dist):
			ax.annotate(row.name, list(row.values),
				xytext=(5,-3), 
				textcoords='offset points',
				size=10, 
				color='darkslategrey')

## create plot of 6 PC combinations
# PC1 vs PC2, PC3 vs PC4 etc.
# arguments are: 
# - pd.DataFrame with PCA transformed gene expression 
# - annotation pd.DataFrame
# - pca sklearn.decomposition object
# - settings object
def plot_2d_pca_multiplot(transformed_expression, annotation, pca, settings):
	fig, ax = plt.subplots(2,3, figsize=(15,10))
	markers = list(annotation["shape"].unique())
	for pc in range(0,12,2):
		for m in markers:
			cells_with_this_shape = annotation["shape"]==m
			ann = annotation.loc[cells_with_this_shape]
			#import pdb; pdb.set_trace()
			transformed_expression.loc[cells_with_this_shape].plot.scatter(
				x=pc+1,
				y=pc+2,
				ax=ax[pc/6][(pc/2)%3],
				s=ann["size"].values,
				c=ann["color"].values,
				legend=True,
				alpha=0.8,
				#edgecolor="black",
				marker = shape_plotly2matplotlib(m)
			)
		
		explained_variance1 = "{0:.2f}".format(pca.explained_variance_ratio_[pc]*100)+"%"
		explained_variance2 = "{0:.2f}".format(pca.explained_variance_ratio_[pc+1]*100)+"%"
		ax[pc/6][(pc/2)%3].set_xlabel("PCA "+str(pc+1)+" ["+explained_variance1+" of variance]")
		ax[pc/6][(pc/2)%3].set_ylabel("PCA "+str(pc+2)+" ["+explained_variance2+" of variance]")
		ax[pc/6][(pc/2)%3].set_aspect("equal", adjustable="box")
	plt.tight_layout()
	plt.subplots_adjust(hspace=0.15, wspace=0.15, left=0.05, bottom=0.05)
	plt.savefig(settings.result_filename+"-pca-multiplot.png", dpi=200)
	plt.show()

## plot cells of defined pair of PCs
# arguments are:
# - pd.DataFrame with PCA transformed gene expression 
# - annotation pd.DataFrame
# - pca sklearn.decomposition object
# - settings object
def plot_2d_pca_single_plot(transformed_expression, annotation, pca, settings, filename=None):
	fig,ax = plt.subplots(figsize=(5,5))
	markers = list(annotation["shape"].unique())
	for m in markers:
		cells_with_this_shape = annotation["shape"]==m
		ann = annotation.loc[cells_with_this_shape]
		transformed_expression.loc[cells_with_this_shape].plot.scatter(
			x=settings.pcs[0],
			y=settings.pcs[1],
			ax=ax,
			s=ann["size"].values,
			c=ann["color"].values,
			legend=True,
			alpha=0.8,
			edgecolor=ann["outline-color"].values,
			marker = shape_plotly2matplotlib(m)
		)
	for cell in transformed_expression.index:
		row = transformed_expression.loc[cell,[int(settings.pcs[0]),int(settings.pcs[1])]]
		df  = transformed_expression.loc[ :  ,[int(settings.pcs[0]),int(settings.pcs[1])]]
		annotate_df(row, df, 8.0, ax)
	
	#ax.set_xlim([-100,100])
	#ax.set_ylim([-100,100])
	ax.set_aspect("equal", adjustable="box")
	plt.xlabel("PCA "+str(settings.pcs[0]))
	plt.ylabel("PCA "+str(settings.pcs[1]))
	plt.tight_layout()
	plt.subplots_adjust(right=0.94)
	if(filename is None):
		filename = settings.result_filename+"PC-"+str(settings.pcs[0])+"-"+str(settings.pcs[1])+".png"
	plt.savefig(filename, dpi=200)
	plt.show()

## record centroids
# arguments are:
# - clusters
# - comb 
# - settings object
def record_trace(clusters, comb, settings, centroids=None):

	test_centroids = centroids #testthis
	
	used_centroids = test_centroids.transpose().iloc[:,[i - 1 for i in settings.pcs]]
	used_centroids.columns = ["x","y","z"]
	used_centroids["color"] = "black"
	used_centroids["shape"] = "shape"
	colors = []
	for i,c in enumerate(clusters):
		#~ ipdb.set_trace()
		colors += [comb.loc[c[1][0],"color"]]
	
	trace = dict(
		name = "centroids",
		x=used_centroids["x"],
		y=used_centroids["y"],
		z=used_centroids["z"],
		type = "scatter3d",
		mode = 'lines+markers',
		line = dict(
			width = 6,
			color = "black", 
			shape = "spline"
		),
		marker = dict(
			size=[15],
			color=used_centroids["color"],
			symbol=["x"]*used_centroids.shape[0], #c[1]["shape"],
			line=dict(width=1) )
		)
	return(trace)

def shape_matplotlib2plotly(s):
		if(s=="o"):
			return "circle"
		elif(s=="s"):
			return "square"
		elif(s=="^"):
			return "triangle-up"
		elif(s=="v"):
			return "triangle-down"
		elif(s=="<"):
			return "triangle-left"
		elif(s==">"):
			return "triangle-right"
		elif(s=="h"):
			return "hexagon"
		elif(s=="*"):
			return "star"
		else:
			return "circle"
			
def shape_plotly2matplotlib(s):
		if(s=="circle"):
			return "o"
		elif(s=="square"):
			return "s"
		elif(s=="triangle-up"):
			return "^"
		elif(s=="triangle-down"):
			return "v"
		elif(s=="triangle-left"):
			return "<"
		elif(s=="triangle-right"):
			return ">"
		elif(s=="hexagon"):
			return "h"
		elif(s=="star"):
			return "*"
		else:
			return "o"

## create 3d PCA plot using plotly library
# arguments are:
# - pd.DataFrame with PCA transformed gene expression 
# - annotation pd.DataFrame
# - settings object
def plot_3d_pca(transformed_expression, annotation, settings, expression_table=None, clusters=None, centroids=None, bin_col_dict=None, height = 1080, width = 1600, genes=None, DEBUG=False):
	used_pcs = transformed_expression[ [settings.pcs[0], settings.pcs[1], settings.pcs[2]]]
	max_range = (used_pcs.max() - used_pcs.min()).max()
	#print(used_pcs.max())
	#print(used_pcs.min())
	layout = dict(
		width=width,
		height=height,
		autosize=False,
		#title='Test',
		scene=dict(
			xaxis=dict(
				# ~ range = [-settings.plot_dim[0], settings.plot_dim[0]],
				range = [used_pcs.iloc[:,0].min()-1, used_pcs.iloc[:,0].min()+max_range+1],
				title="PC "+str(settings.pcs[0]),
				gridcolor='rgb(0, 0, 0)',
				zerolinecolor='rgb(255, 0, 0)',
				showbackground=True,
				backgroundcolor='#bababa'
			),
			yaxis=dict(
				range = [-settings.plot_dim[1], settings.plot_dim[1]],
				# ~ range = [used_pcs.iloc[:,1].min()-1, used_pcs.iloc[:,1].min()+max_range+1],
				title="PC "+str(settings.pcs[1]),
				gridcolor='rgb(0, 0, 0)',
				zerolinecolor='rgb(255, 0, 0)',
				showbackground=True,
				backgroundcolor='#bababa'
			),
			zaxis=dict(
				range = [-settings.plot_dim[2], settings.plot_dim[2]],
				# ~ range = [used_pcs.iloc[:,2].min()-1, used_pcs.iloc[:,2].min()+max_range+1],
				title="PC "+str(settings.pcs[2]),
				gridcolor='rgb(0, 0, 0)',
				zerolinecolor='rgb(255, 0, 0)',
				showbackground=True,
				backgroundcolor='#bababa'
			),
			aspectmode = 'manual'
		),
	)
	comb = pd.concat([transformed_expression, annotation], axis=1)
	#comb["name"] = comb["shape"]
	if "name" not in comb.columns:
		comb["name"] = comb["color"]+"_"+comb["shape"]
	traces = comb["name"].unique()
	# allow coloring cells by quantile expression of supplied genes
	if (genes is not None):
		markers = list(annotation["shape"].unique())
		mg = mygene.MyGeneInfo()
		# ~ for i in enumerate(genes): 
		gene_info = mg.querymany(genes, scopes='symbol', fields='ensembl.transcript')[0]
		if len(gene_info['ensembl']) > 1:
			trx = gene_info['ensembl'][0]['transcript']
		else:
			trx = gene_info['ensembl']['transcript']
		#test = all(elem in test for elem in expression_table.columns)
		if type(trx) == list:
			sum_trx = expression_table.loc[:,trx].sum(axis=1)
			# catch pandas omission of missing list values in loc; incompatible between pandas versions
			#~ sub_trx = expression_table.reindex(trx, axis = 1)
			#~ sum_trx = sub_trx.loc[:,trx].sum(axis=1)
		else:
			sum_trx = expression_table.loc[:,trx]
		
		# color by quantile
		trx_df = sum_trx.to_frame('trx_count')
		number_of_bins = len(bin_col_dict.keys())
		bin_labels=bin_col_dict.keys()
		
		trx_df['bin'] = pd.cut(trx_df.trx_count, number_of_bins, labels=bin_labels)
		trx_df['color'] = trx_df['bin'].map(bin_col_dict)
		comb['color'] = trx_df['color']
		comb['name'] = trx_df['bin']
		traces = comb["name"].unique().sort_values()

	data = []
	for t in traces:

		trace = dict(
			name=t,
			text = comb.loc[comb["name"]==t,:].index,# + " "+ transformed_expression["day"], #+ "\n" + transformed_expression["branch"],
			x = comb.loc[comb["name"]==t,settings.pcs[0]],
			y = comb.loc[comb["name"]==t,settings.pcs[1]],
			z = comb.loc[comb["name"]==t,settings.pcs[2]],
			type = "scatter3d",
			mode = 'markers',
			marker = dict(
				size=comb.loc[comb["name"]==t,"size"].values,
				# ~ color=trx_df.color,
				color=comb.loc[comb["name"]==t,"color"].values,
				symbol=comb.loc[comb["name"]==t,"shape"].apply(shape_matplotlib2plotly).values,
				line=dict(width=1) )
			)
		data.append( trace )
	
	if(clusters != None):
		
		centroids = get_cluster_centroids(transformed_expression, clusters)
		trace = record_trace(clusters, comb, settings, centroids)
		data.append(trace)
	if(DEBUG):
		IPython.embed()
		
		
	#print(annotation["shape"].apply(shape_matplotlib2plotly))
	
	# check if start/end centroids extend off the plot dimensions and reassign dimensions if so
	def nested_set(dic, keys, value):                     
		for key in keys[:-1]:
			dic = dic.setdefault(key, {})
		dic[keys[-1]] = value     
	
	axis_tuples = list(zip(["x", "y", "z"],["xaxis", "yaxis", "zaxis"]))
	for i,c in axis_tuples:
		if (trace[i].min() < layout['scene'][c]['range'][0]):
			nested_set(layout, ['scene', c, 'range'], [trace[i].min(),layout['scene'][c]['range'][1]])
		if (trace[i].max() > layout['scene'][c]['range'][1]):
			nested_set(layout, ['scene', c, 'range'], [layout['scene'][c]['range'][0],trace[i].max()])
			
	# ~ fig = dict(data=data, layout=layout)
	if (genes is not None):
		fig = go.Figure(data=data, layout=layout)
		url = plotly.offline.plot(fig, filename=settings.result_filename+"_"+genes, validate=False, auto_open=False)
	else:
		fig = dict(data=data, layout=layout)
		url = plotly.offline.plot(fig, filename=settings.result_filename, validate=False, auto_open=False)


	
	print(url)
	return(fig)

## set pca dimension in nested dictionary object
#  takes as input:
#  - dictionary object
#  - list of nested keys ['scene', 'xaxis', 'range']
#  - new value to be set layout['scene']['xaxis']['range'][0],0]	
def nested_set(dic, keys, value):                     
		for key in keys[:-1]:
			dic = dic.setdefault(key, {})
		dic[keys[-1]] = value 

## take linkage (from scipy.linkage) and generate n clusters
#  takes as input:
#  - linkage
#  - number of clusters to generate
#  - labels of nodes (linkage contains only ordinal number of each node)
def get_cluster_labels(linkage, n_clusters, labels):
	n = labels.shape[0]
	clusters = {}
	for i in range(0,n-n_clusters):
		if(linkage[i,0]<n):
			a=set([labels[int(linkage[i,0])]])
		else:
			a=clusters[int(linkage[i,0])]
			del(clusters[int(linkage[i,0])])
		if(linkage[i,1]<n):
			b=set([labels[int(linkage[i,1])]])
		else:
			b=clusters[int(linkage[i,1])]
			del(clusters[int(linkage[i,1])])
		clusters[n+i] = a.union(b)
	return clusters

## change colors of cells according to the clusters provided
#  clusters is a dictionary cluster_number => set of cells
#  annotation is a pandas dataframe as generated by read_expression
#  colors are the colors to assign
def change_annotation_colors_to_clusters(clusters, annotation, colors):
	for i,c in enumerate(clusters.values()):
		annotation.loc[c, "color"] = colors[i]
	print colors

## plot hierarchical clustering for all methods of linkage
# arguments are:
# - pd.DataFrame with PCA transformed gene expression 
# - annotation pd.DataFrame
# - settings object
# - filename for output picture
def plot_all_hierarchical_clusterings(transformed_expression, annotation, settings, clusters=None):
	
	link_color = {}
	def link_color_func(node):
		return link_color[node]
	
	scipy_linkage_methods = ["complete", "average", "single", "centroid", "median", "ward"]
	# ~ scipy_linkage_methods = ["ward"]
	# plot clusterings on one magor figure
	#fig,ax = plt.subplots(nrows=2, ncols=3, figsize=(50, 30))
	i=0
	pp = PdfPages(settings.result_filename+"-clustering.pdf")
	for method in scipy_linkage_methods:
		plot_hierarchical_clustering(transformed_expression, annotation, method=method, sett=settings, clusters=clusters)
		i += 1
		pp.savefig()
	pp.close()
	print("save results to:\n"+settings.result_filename+"-clustering.png"+"\n"+settings.result_filename+"-clustering.pdf")
	plt.savefig(settings.result_filename+"-clustering.png", dpi=200)


## plot hierarchical clustering with a given linkage method
#  used when assigning centroids based on clustering
#  arguments are:
# - pd.DataFrame with PCA transformed gene expression
# - method of linkages (ward, complete, average, etc.)
def plot_hierarchical_clustering(transformed_expression, annotation, method, color_scheme="static", sett=settings, clusters=None):
	# color links on the basis of connection to same-group neighbor. 
	# If neighbors in same group, color identically. If neighbors in different groups, color gray.
	def colorize_links(linkage):
		l_color = {}
		n = transformed_expression.shape[0]
		for i in range(0,n):
			l_color[i] = annotation.iloc[i,]["color"]
		#print l_color
		
		# ~ test = sch.fcluster(linkage, sett.num_clusters, criterion='maxclust')
		# ~ cluster_col_dict = dict(zip(range(1,sett.num_clusters+1),["blue", "green", "red", "yellow"]))
		# ~ test2 = pd.DataFrame(test)
		# ~ test2['color'] = test2[0].map(cluster_col_dict)
		for i in range(0,linkage.shape[0]):
			clust1 = int(linkage[i,0])
			clust2 = int(linkage[i,1])
			#print clust1, clust2
			if(l_color[clust1] == l_color[clust2]):
				l_color[n+i] = l_color[clust1]
			elif (color_scheme == "dynamic"):
				l_color[n+i] = "gray"
			else:
				l_color[n+i] = l_color[clust2]
		#print l_color
		return l_color
	
	def link_color_func(node):
		return link_color[node]
			
	linkage = sc.cluster.hierarchy.linkage(transformed_expression, method=method)
	
	
	clusters_without_time = get_cluster_labels(linkage, sett.num_clusters, transformed_expression.index)
	cluster_colors = ["blue", "red", "orange", "purple", "green", "brown", "black", "gray", "lawngreen", "magenta", "lightpink", "indigo", "lightblue", "lightgoldenrod1", "mediumpurple2"]
	change_annotation_colors_to_clusters(clusters_without_time, annotation, cluster_colors)
	clusters = []
	
	for i in range(0,sett.num_clusters):
		time = i+1
		clusters.append( (time,annotation.loc[annotation["color"]==cluster_colors[i]].index) )
	clusters.sort(key=lambda by_first: by_first[0])
	
	link_color = colorize_links(linkage)
	fig,ax = plt.subplots(figsize=(32,10))
	dendro  = sc.cluster.hierarchy.dendrogram(
		linkage,
		#ax=ax[i/3,i%3],
		ax=ax,
		labels = transformed_expression.index,
		link_color_func = link_color_func,
		#color_threshold = 0,
		#above_threshold_color = "black",
		count_sort = "ascending") #, title=method
	#ax[i/3,i%3].set_title(method)
	ax.set_title(method)
	#tick_labels = ax[i/3,i%3].get_xmajorticklabels()
	tick_labels = ax.get_xmajorticklabels()
	for lbl in tick_labels:
		lbl.set_color(annotation.loc[lbl.get_text()]["color"])
	#return(dendro)

## rotate transformed expression matrix by defined angle
#  used internally in order to define pseudotime
#  arguments are:
# - pd.DataFrame with PCA transformed gene expression 
# - x,y = PCs to rotate
# - angle in degrees
# returns:
# pdDataFrame with values in columns x,y rotated by angle
def rotate_expression(transformed_expression,x,y,angle):
	theta = math.radians(angle)
	ret = transformed_expression.copy()
	ret[x] = transformed_expression[x]*math.cos(theta) - transformed_expression[y]*math.sin(theta)
	ret[y] = transformed_expression[x]*math.sin(theta) + transformed_expression[y]*math.cos(theta)
	return ret

## function 
# - finds pair of 2 PCs that are most correlated with time labels (as defined by "day" column in annotation table) using spearman correlation
# - finds rotation of this PCs so X axis has best correlation with time
# 
# returns: pseudotime for each cell, defined as linear combination of PCs, having best time correlation
# 
# arguments are:
# - pd.DataFrame with PCA transformed gene expression 
# - annotation pd.DataFrame
# - pca sklearn.decomposition object
# - settings object
def find_pseudotime(transformed_expression, annotation, pca, settings, user_pcs=None):
	n_pca = len(transformed_expression.columns)
	transformed_expression["day"] = annotation["day"]
	transformed_expression_without_superimposed = transformed_expression.loc[annotation[annotation["superimpose-for-spearman"]==False].index]
	print "Finding best rotation for Spearman correlation. Shape of used table:",transformed_expression_without_superimposed.shape
	spearman = transformed_expression_without_superimposed.corr(method="spearman").loc["day",range(1,n_pca+1)].abs().sort_values(ascending=False)
	#plot_spearman correlations and explained variation
	spearman_filename = settings.result_filename.replace(".png", "_spearman.png")
	width=0.2
	
	fig,ax = plt.subplots(figsize=(8,5))
	
	ax2= ax.twinx()
	
	spearman.plot.bar(ax=ax, width=width, position=1, color="blue")
	pd.Series(pca.explained_variance_ratio_, index=range(1,n_pca+1)).loc[spearman.index].plot.bar(ax=ax2, width=width, position=0, color="red")
	
	ax.set_xlabel("PC component")
	ax.set_ylabel("Spearman correlation\nto days [blue]")
	ax2.set_ylabel("% variance explained [red]")
	
	
	# plot difference between shRBKD and shCtrl cells if both present
	# calculate distance between RBKD and Ctrl cells
	shCtrl_cells = annotation.loc[annotation['treatment'].str.contains("shCtrl")].index
	vals_shCtrl = transformed_expression.loc[shCtrl_cells,:]

	if len(shCtrl_cells) > 0:
		shRBKD_cells = annotation.drop(shCtrl_cells).index
		vals_shRBKD = transformed_expression.loc[shRBKD_cells,:]
		t_val = []
		for i in transformed_expression.iloc[:, :-1]:
			test = (abs(vals_shCtrl.loc[:,i].mean()-vals_shRBKD.loc[:,i].mean()))/np.std(vals_shCtrl.loc[:,i].append(vals_shRBKD.loc[:,i]))
			t_val.append(test) 
		
		#~ # define axis 3
		ax3= ax.twinx()
		#~ # assign location for axis 3
		ax3.spines['right'].set_position(('outward', 60))
		#~ # plot axis 3 green bars
		pd.Series(t_val, index=range(1,n_pca+1)).loc[spearman.index].plot.bar(ax=ax3, width=width, position=2, color="green")
		#~ # assign label for axis 3
		ax3.set_ylabel("distance between RBKD and Ctrl [green]")
	plt.tight_layout()
	low,high = plt.xlim()
	plt.xlim(low-0.5, high)
	plt.savefig(spearman_filename, dpi=200)
	#~ ipdb.set_trace()
	if user_pcs:
		settings.pcs = user_pcs
	else:
		settings.pcs = spearman.iloc[0:2].index
	
	# find best rotation
	best_angle = 0
	best_spearman = 0
	for a in range(0,360):
		te = rotate_expression(transformed_expression_without_superimposed, settings.pcs[0], settings.pcs[1], a)
		spearman = te.corr(method="spearman").loc["day",int(settings.pcs[0])]
		#print "Trying angle: ",a," spearman: ",spearman
		if(spearman > best_spearman):
			best_angle = a
			best_spearman = spearman
		
	del(transformed_expression["day"])
	print settings.pcs
	print "Best rotation: ",best_angle
	
	rotated_expression = rotate_expression(transformed_expression, int(settings.pcs[0]), int(settings.pcs[1]), best_angle)
	# plot original PC plot
	plot_2d_pca_single_plot(transformed_expression, annotation, pca, settings, filename = settings.result_filename+"-original")
	# plot rotated PC plot
	plot_2d_pca_single_plot(rotated_expression, annotation, pca, settings, filename = settings.result_filename+"-rotated")
	pt = rotated_expression[int(settings.pcs[0])]
	pt.name = "pseudotime"
	# normalize to <0;1>
	pt = (pt-pt.min())/(pt.max()-pt.min())
	print("printing correlation plot to "+spearman_filename+".png")
	return pt

## function reads list of integers and/or ranges
#  and converts it to list of integers
#  "1,3-5,8" becomes [1,3,4,5,8]
def list_from_ranges(s):
	l = []
	for i in s.split(","):
		j = i.split("-")
		if(len(j)==1):
			l.append(int(i))
		elif(len(j)==2):
			x,y = int(j[0]), int(j[1])
			for k in range(x,y+1):
				l.append(k)
	return l

## function calculates pseudotime for each cell based on defined times for clusters
#  pseudotime for each cell is calculated as weighted average of times assigned to clusters
#  weight for each cluster is an inverse square distance of the cell to the cluster
#  w = 1/(dist^2)
def calculate_pseudotime_using_cluster_times(PC_expression, annotation, clusters, settings):
	
	palette_size = int(raw_input("What palette size would you like to use (how many colors)? "))
	calculate_on = list_from_ranges(raw_input("Which PCs would you like to use for calculating pseudotime? [type comma separated list, list can also include ranges 1-5,8] "))
	used_PC_expression = PC_expression[calculate_on]
	centroids = get_cluster_centroids(used_PC_expression, clusters)
	sq_distances = pd.DataFrame(index=used_PC_expression.index, columns=[])
	weights = pd.DataFrame(index=used_PC_expression.index, columns=[])
	test = pd.DataFrame(index=used_PC_expression.index, columns=[])
	for i,c in enumerate(centroids):
		sq_distances[i] = ((used_PC_expression-centroids[c])**2).sum(axis=1)**0.5
		weights[i] = 1/sq_distances[i]
		test += sq_distances[i]
	
	pseudotime_clusters = [(clusters[0][0]-1,clusters[0][1])]
	pseudotime_clusters.extend(clusters)
	pseudotime_clusters.append(tuple((clusters[-1][0]+1, clusters[-1][1])))
	
	# fine tune weighting
	#~ c_weights, c_clusts =  map(list, zip(*pseudotime_clusters))
	#~ new_weights = [x for x in c_weights]
	#~ pseudotime_clusters = zip(new_weights, c_clusts)
	
	# exclude first and last "real centroids"
	#~ cols = [1,-2]
	#~ weights.drop(weights.columns[cols], axis=1, inplace = True)
	#~ weights.columns = range(0,weights.shape[1])

	pseudotime = pd.Series(0, index=used_PC_expression.index)
	for w in weights:
		print w
		pseudo_part = (pseudotime_clusters[w][0]+1)*weights[w]
		pseudotime += pseudo_part
	pseudotime /= weights.sum(axis=1)
	# normalize
	pseudotime -= pseudotime.min()
	pseudotime /= pseudotime.max()
	pal = sns.cubehelix_palette(palette_size+1, start=2, rot=0, dark=0, light=0.85)
	pal = [(int(i[0]*256),int(i[1]*256),int(i[2]*256)) for i in pal]
	color_indices = map(int,pseudotime*palette_size)
	annotation["color"] = [RGBToHTMLColor(pal[i]) for i in color_indices]
	return pseudotime, centroids


## fits polynomial curve on gene expression data, and returns value of this curve in equal intervals over pseudotime
# arguments:
# - pd.DataFrame with gene expression 
# - pd.Series with pseudotime coordinates for each cell
# - Ensamble transcript ID
# - degree of the polynomial to fit
# - number of samples to return (they will be equally spaced over pseudotime)
def interpolate_gene_over_pseudotime(exp, pseudotime, transcript_id, weights=None, degree=3, n_samples=20):
	#expr_over_ptime = pd.DataFrame(pseudotime)
	#expr_over_ptime["expression"] = exp[transcript_id]
	curve_coeff, res, _, _, _ = np.polyfit(x = pseudotime, y=exp[transcript_id], deg = degree, full=True, w=weights)
	curve_func  = np.poly1d(curve_coeff)
	samples = np.linspace(0,1,n_samples)
	fitted_curve = pd.DataFrame([(time, curve_func(time)) for time in samples], columns = ["pseudotime", "expression"])
	return fitted_curve,res

## plots gene expression over pseudotime
# arguments are:
# - pd.DataFrame with gene expression 
# - pd.Series with pseudotime coordinates for each cell
# - Ensamble transcript ID
def plot_gene_with_pseudotime(exp, pseudotime, transcript_id, annotation, filename=None, ax=None, plot_id=None, ctrl_pseudotime=None):
	expr_over_ptime = pd.DataFrame(pseudotime)
	#~ IPython.embed()
	expr_over_ptime["expression"] = exp.loc[pseudotime.index, transcript_id]
	if ctrl_pseudotime is not None:
		ctrl_over_ptime = pd.DataFrame(ctrl_pseudotime)
		ctrl_over_ptime["expression"] = exp.loc[ctrl_pseudotime.index, transcript_id]
	
	# translast colors by day (in ctrl cells)
	day_list = list(annotation.day.unique())
	color_list = list(annotation.color.unique())
	if 'grey' in color_list:
		color_list.remove('grey')
	color_by_day = dict(zip(day_list,color_list))
	def day_to_color(row, color_day_dict):
		return(color_day_dict[row['day']])

	if plot_id == "exp":
		expr = annotation.loc[(~annotation.treatment.str.contains('shCtrl')),:]
		RBKD_over_ptime = expr_over_ptime[expr_over_ptime.index.isin(expr.index)]

		expr_ann = annotation.loc[RBKD_over_ptime.index, :] 
		# ~ ax = RBKD_over_ptime.plot.scatter(x="pseudotime", y="expression", c=expr_ann["color"], ax=ax)
		ax = RBKD_over_ptime.plot.scatter(x="pseudotime", y="expression", c=expr_ann["color"], ax=ax)
		lowess = sm.nonparametric.lowess
		z = lowess(RBKD_over_ptime["expression"], pseudotime[pseudotime.index.isin(RBKD_over_ptime.index)])
		pd.DataFrame(z, columns=["pseudotime","local regression"]).plot.line(x="pseudotime", y="local regression", c="gray", style="--", ax=ax)
	elif plot_id == "Ctrl_wo_RBKD":
		shctrl = annotation.loc[(annotation.treatment.str.contains('shCtrl')),:]
		ctrl_over_ptime = expr_over_ptime[expr_over_ptime.index.isin(shctrl.index)]

		ctrl_ann = annotation.loc[ctrl_over_ptime.index, :] 
		
		# convert ctrl cells from gray to colored for plotting 
		translate_colors = ctrl_ann.apply(day_to_color, args=(color_by_day,), axis=1)
		ax = ctrl_over_ptime.plot.scatter(x="pseudotime", y="expression", c=translate_colors, ax=ax)
		lowess = sm.nonparametric.lowess
		#~ IPython.embed()
		z = lowess(ctrl_over_ptime["expression"], pseudotime[pseudotime.index.isin(ctrl_over_ptime.index)])
		subplt = pd.DataFrame(z, columns=["pseudotime","local regression"]).plot.line(x="pseudotime", y="local regression", c="gray", style="--", ax=ax)
		# hardcode x-axis so that control is directly comparable to RBKD
		subplt.set_xlim(0.0,1.0)
	elif plot_id == "Ctrl_alone":
		shctrl = annotation.loc[(annotation['treatment']=="shCtrl"), :]
		ctrl_over_ptime = ctrl_over_ptime[ctrl_over_ptime.index.isin(shctrl.index)]

		ctrl_ann = shctrl.loc[ctrl_over_ptime.index, :] 
		
		# convert ctrl cells from gray to colored for plotting 
		translate_colors = ctrl_ann.apply(day_to_color, args=(color_by_day,), axis=1)
		
		ax = ctrl_over_ptime.plot.scatter(x="pseudotime", y="expression", c=translate_colors, ax=ax)

		lowess = sm.nonparametric.lowess
		z = lowess(ctrl_over_ptime["expression"], ctrl_pseudotime)
		pd.DataFrame(z, columns=["pseudotime","local regression"]).plot.line(x="pseudotime", y="local regression", c="gray", style="--", ax=ax)

	
	# ~ ax.legend_.remove()
	#plt.tight_layout()
	if(filename==None):
		#plt.show()
		pass
	else:
		plt.savefig(filename)
		plt.close('all')

## read pseudotime from tab delimited csv file
def read_pseudotime_from_file(filename):
	return pd.read_csv(filename, sep="\t", index_col=0, names=["pseudotime"])["pseudotime"]

## look up gene, transcript dict from mygene
def get_gene_transcript_dic(expression_table):
	transcripts = expression_table.columns.copy()
			
	mg = mygene.MyGeneInfo()

	gene_info = mg.querymany(transcripts, scopes='ensembl.transcript', fields='symbol', returnall=False)
	
	# remove transcripts not found in gene lookup
	gene_info[:] = [d for d in gene_info if d.get('notfound') != True]
	
	# define contiguous list of transcripts
	transcripts = [query for i, query in enumerate(d['query'] for d in gene_info)]
	# define contiguous list of genes (matching transcripts)
	genes = [symbol for i, symbol in enumerate(d['symbol'] for d in gene_info)]
	
	dic={}
	for x,y in zip(transcripts, genes):
		dic.setdefault(y,[]).append(x)	
	
	return(dic)

## look up gene, transcript dict from mygene
def trx_to_gene_exp_table(expression_table, gene_trx_dic):
	gene_exp = []
	for i,gene in enumerate(gene_trx_dic):
	#~ for i,gene in enumerate(dic):
		if i%1000 == 0:
			#~ print gene
			print "Genes processed:",i
		gene_col = expression_table.loc[:, gene_trx_dic[gene]].sum(axis=1)
		gene_col.columns = [gene]
		#generate gene-level exprescsion table
		gene_exp.append(gene_col)
		#~ IPython.embed()
	
	gene_exp = pd.concat(gene_exp, axis = 1)
	gene_exp.columns = gene_trx_dic.keys()

	return(gene_exp)

## returns spearman correlation of each gene in expression matrix with pseudotime
# arguments are:
# - exp = pd.DataFrame with gene expression 
# - pseudotime = pd.Series with pseudotime coordinates for each cell
# - [optional] correlation_threshold = returns only genes with absolute value of correlation >= threshold
def get_correlation_with_pseudotime(pseudotime, exp, annotation, gene_trx_dic, cell_set_flag=None, feature = "gene", correlation_threshold = 0, method = "spearman"):
	
	def return_subset_correlation(subset_index, feature):
			subset_index = pseudotime.index[pseudotime.index.isin(subset_index)]
			transcripts = exp.columns.copy()
			subsetc = exp.loc[subset_index]
			subsetc["pseudotime"] = pseudotime[subset_index]
			if feature == "gene":
				spearman = pd.DataFrame(0, index=gene_trx_dic.keys(), columns=["corr"])
				# correlation by gene
				for i,gene in enumerate(gene_trx_dic):
				#~ for i,gene in enumerate(gene_trx_dic):
					if i%1000 == 0:
						#~ print gene
						print "Genes processed:",i
					gene_col = subsetc.loc[:, gene_trx_dic[gene]].sum(axis=1)
					gene_col.columns = [gene]
					corr = pd.concat([gene_col, subsetc.loc[:,"pseudotime"]], axis=1)
					corr.columns = [gene, 'pseudotime']
					corr = corr.loc[ : , [gene,"pseudotime"]].corr(method=method).iloc[0,1]
					if corr != corr: # if NaN (no data to calculate on)
						corr = 0 # then correlation is zero
					spearman.loc[gene,"corr"] = corr
			elif feature == "transcript":
				spearman = pd.DataFrame(0, index=transcripts, columns=["corr"])
				# correlation by transcript
				for i,transcript in enumerate(transcripts):
					if i%1000 == 0:
						print "Transcripts processed:",i
					corr = subsetc.loc[ : , [transcript,"pseudotime"]].corr(method=method).iloc[0,1]
					if corr != corr: # if NaN (no data to calculate on)
						corr = 0 # then correlation is zero
					spearman.loc[transcript,"corr"] = corr
			return(spearman)
	
	
	if cell_set_flag == "ctrl":
		spearman = return_subset_correlation(pseudotime.index, feature)

	elif cell_set_flag == "exp":
		spearman = return_subset_correlation(pseudotime.index, feature)

	else:
		exp_index = annotation.loc[annotation["treatment"]!="shCtrl"].index
		shctrl_index = annotation.loc[annotation["treatment"]=="shCtrl"].index
		# ~ print(cell_set_flag)
		# ~ print(shctrl_index)
		if shctrl_index.empty:
			subset_indices = [exp_index]
			cell_set_flags = ["exp"]
			# check if map is returning spearman correlation and gene_expression_table
			spearman = [return_subset_correlation(x, feature) for x in subset_indices]
			#~ spearman = map(return_subset_correlation, subset_indices, feature)
			spearman = pd.concat(spearman, axis=1)
			spearman.columns = cell_set_flags
		else:
			subset_indices = [exp_index, shctrl_index]
			cell_set_flags = ["RBKD", "shCtrl"]
			# check if map is returning spearman correlation and gene_expression_table
			spearman = [return_subset_correlation(x, feature) for x in subset_indices]
			#~ spearman = map(return_subset_correlation, subset_indices, feature)
			spearman = pd.concat(spearman, axis=1)
			spearman.columns = cell_set_flags
			#~ IPython.embed()

	return spearman
	

def plot_3d_pca_colored_by_clustering(PC_expression, annotation, pca, settings):
	
	link_color = {}
	def link_color_func(node):
		return link_color[node]
	
	def colorize_links(linkage):
		l_color = {}
		n = PC_expression.shape[0]
		for i in range(0,n):
			l_color[i] = annotation.iloc[i,]["color"]
		#print l_color
		for i in range(0,linkage.shape[0]):
			clust1 = int(linkage[i,0])
			clust2 = int(linkage[i,1])
			#print clust1, clust2
			if(l_color[clust1] == l_color[clust2]):
				l_color[n+i] = l_color[clust1]
			else:
				l_color[n+i] = "gray"
		#print l_color
		return l_color
	
	palette = ["red","blue","green","orange","purple","pink","#f47442"]
	linkage = sc.cluster.hierarchy.linkage(PC_expression[settings.pcs], method=settings.clustering_method)
	clusters = get_cluster_labels(linkage, settings.n_clusters, PC_expression.index)
	csv_filename = settings.result_filename+"PCs-"+"-".join(map(str,settings.pcs))+"-clustered-"+str(settings.n_clusters)+".csv"
	with open(csv_filename,"w") as f:
		for i,c in enumerate(clusters):
			f.write("cluster_"+str(i)+"\t")
			f.write("\t".join(clusters[c]))
			f.write("\n")
			annotation.loc[clusters[c],"color"] = palette[i]
	
	link_color = colorize_links(linkage)
	
	plot_3d_pca(PC_expression, annotation, settings)
	dendro  = sc.cluster.hierarchy.dendrogram(linkage,labels = PC_expression.index,count_sort = "ascending", link_color_func = link_color_func)
	plt.show()

## takes annotation dataframe and returns list of day clusters
#  as list of toupples (day, index_of_cells)
def time_clusters_from_annotations(annotation):
	days = annotation["day"].sort_values().unique()
	clusters = []
	for d in days:
		clusters.append((d,annotation[annotation["day"]==d].index))
	return clusters

## takes PCA transformed expression and list of clusters [(time, index_of_cells), ...]
#  and returns centroid for each cluster
def get_cluster_centroids(PC_expression, clusters):
	centroids = []
	for cl in clusters:
		centroids.append(PC_expression.loc[cl[1],:].mean())
	# append "first cell" and "last cell" to centroids to 
	# based on extrapolation of trace from centroids 0 to 1 and n-1 to n respectively
	half_trace_seg_1 = (centroids[0] - centroids[1])/2
	half_trace_seg_n = (centroids[-1] - centroids[-2])/2
	first_cell = centroids[0] + half_trace_seg_1
	last_cell = centroids[-1] + half_trace_seg_n
	centroids.insert(0, first_cell)
	centroids.append(last_cell)
	
	centroids = pd.concat(centroids, axis=1)
	return centroids

def plot_heatmap(expression_table, annotation, cell_dendro):
	# get data
	data = np.genfromtxt("http://files.figshare.com/2133304/ExpRawData_E_TABM_84_A_AFFY_44.tab",
						 names=True,usecols=tuple(range(1,30)),dtype=float, delimiter="\t")
	data_array = data.view((np.float, len(data.dtype.names)))
	data_array = data_array.transpose()
	labels = data.dtype.names
	
	expression_table = expression_table.loc[cell_dendro["ivl"]]
	data_array0 = expression_table.as_matrix()
	data_array0 = data_array0.transpose()
	labels0 = expression_table.index
	
	
	#~ # Initialize figure by creating side dendrogram
	#~ figure = FF.create_dendrogram(data_array0, orientation='right')
	#~ for i in range(len(figure['data'])):
		#~ figure['data'][i]['xaxis'] = 'x2'
	
	# Create Upper Dendrogram
	figure = FF.create_dendrogram(data_array, orientation='top', labels=labels)
	for i in range(len(figure['data'])):
		figure['data'][i]['yaxis'] = 'y2'
	
	#~ dendro_top = FF.create_dendrogram(data_array, orientation='top', labels=labels)
	#~ for i in range(len(dendro_top['data'])):      
		#~ dendro_top['data'][i]['yaxis'] = 'y2'
	
	# Add Top Dendrogram Data to Figure
	#~ figure['data'].extend(dendro_top['data'])

	#convert scipy dend to plotly
	
	cell_dendro_k, cell_dendro_v = cell_dendro.keys(), cell_dendro.values()
	
	cell_ids = cell_dendro_v[0]
	
	ptl_dend = dict(zip(cell_ids, zip(cell_dendro_v[1], cell_dendro_v[2], cell_dendro_v[3], cell_dendro_v[4])))
	
	new_dend = {}
	for key,value in ptl_dend.iteritems():
		layout = dict(
			hoverinfo = 'text',
			marker =  dict(color = value[2]),
			mode = 'lines',
			text = None,
			type = 'scatter',
			x = value[3],
			xaxis = 'x2',
			y = value[0],
			yaxis = 'y')
		new_dend[key] = layout

	def without_keys(d, keys):       
		return {x: d[x] for x in d if x not in keys}
		
	new_top = without_keys(figure, "data")
	
	new_top["layout"]["xaxis"]["ticktext"] = ptl_dend.keys()

	new_tickvals = [5.0 + 10*idx for idx,val in enumerate(ptl_dend.keys())]
	
	new_top["layout"]["xaxis"]["tickvals"] = new_tickvals

	new_top["data"] = new_dend

	# Add Side Dendrogram Data to Figure
	#~ figure['data'].extend(new_top)
	
	# Create Heatmap
	dendro_leaves = new_top['layout']['xaxis']['ticktext']
	dendro_leaves = [i for i,v in enumerate(dendro_leaves)]

	heat_data = pdist(data_array0)
	heat_data = squareform(data_dist)
	
	#~ heat_data = heat_data[dendro_leaves,:]
	#~ heat_data = heat_data[:,dendro_leaves]

	IPython.embed()

	heatmap = Data([
		Heatmap(
			z = heat_data,
			x = dendro_leaves,
			y = expression_table.index,
			colorscale = 'YIGnBu'
		)
	])
	
	IPython.embed()
	
	heatmap[0]['x'] = figure['layout']['xaxis']['tickvals']
	#~ heatmap[0]['y'] = dendro_side['layout']['yaxis']['tickvals']

	# Add Heatmap Data to Figure
	figure['data'].extend(Data(heatmap))

	# Edit Layout
	figure['layout'].update({'width':800, 'height':800,
							 'showlegend':False, 'hovermode': 'closest',
							 })
	# Edit xaxis
	figure['layout']['xaxis'].update({'domain': [.15, 1],
									  'mirror': False,
									  'showgrid': False,
									  'showline': False,
									  'zeroline': False,
									  'ticks':""})
	# Edit xaxis2
	figure['layout'].update({'xaxis2': {'domain': [0, .15],
									   'mirror': False,
									   'showgrid': False,
									   'showline': False,
									   'zeroline': False,
									   'showticklabels': False,
									   'ticks':""}})

	# Edit yaxis
	figure['layout']['yaxis'].update({'domain': [0, .85],
									  'mirror': False,
									  'showgrid': False,
									  'showline': False,
									  'zeroline': False,
									  'showticklabels': False,
									  'ticks': ""})
	# Edit yaxis2
	figure['layout'].update({'yaxis2':{'domain':[.825, .975],
									   'mirror': False,
									   'showgrid': False,
									   'showline': False,
									   'zeroline': False,
									   'showticklabels': False,
									   'ticks':""}})
	url = plotly.offline.plot(figure, filename='dendrogram_with_heatmap', validate=False, auto_open=False)
	print(url)

## main function
#  when run separately, program expects following arguments:
# - argv[1] = comma separated file with expression values
# - argv[2] = file with cell sets (see settings.read_cell_sets())
# - argv[3] = file with commands for program (what to do). Usually named plot_settings.csv
def main():
	# get parameters
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
	#print "Running in mode:",sett.run_mode
	
	if(sett.run_mode=="2d-pca-multiplot"):
		plot_2d_pca_multiplot(PC_expression, annotation, pca, sett)
	elif(sett.run_mode=="2d-pca-single"):
		plot_2d_pca_single_plot(PC_expression, annotation, pca, sett)
	elif(sett.run_mode=="3d-pca"):
		plot_3d_pca(PC_expression, annotation, sett)
	elif(sett.run_mode=="hierarchy"):
		plot_all_hierarchical_clusterings(PC_expression, annotation, sett)
	elif(sett.run_mode=="pseudotime"):
		pseudotime = find_pseudotime(subset_PC_expression, subset_annotation, pca, sett) #comment out
		pseudotime.to_csv(sett.result_filename+"_pseudotime.csv", sep="\t") #comment out
		time_clusters_from_annotations(annotation) #formerly get_time_clusters_from_annotations
	elif(sett.run_mode == "3d-pca-colored-by-clustering"):
		plot_3d_pca_colored_by_clustering(PC_expression, annotation, pca, sett)
		
	elif(sett.run_mode == "test"):
		palette_size = 10
		clusters = time_clusters_from_annotations(annotation)
		print(clusters)
		calculate_on = [1,2,3]
		used_PC_expression = PC_expression[calculate_on]
		centroids = get_cluster_centroids(used_PC_expression, clusters)
		print(centroids)
		sq_distances = pd.DataFrame(index=used_PC_expression.index, columns=[])
		weights = pd.DataFrame(index=used_PC_expression.index, columns=[])
		for i,c in enumerate(centroids):
			sq_distances[i] = ((used_PC_expression-c)**2).sum(axis=1)**0.5
			weights[i] = 1/sq_distances[i]
		
		print(weights[0])
		pseudotime = pd.Series(0, index=used_PC_expression.index)
		for w in weights:
			print(clusters[w][0])
			print(w)
			pseudotime += clusters[w][0]*weights[w]
		pseudotime /= weights.sum(axis=1)
		# normalize
		pseudotime -= pseudotime.min()
		pseudotime /= pseudotime.max()
		#pseudotime += 0.1
		print(pseudotime)
		#annotation["size"] = (pseudotime*10)+3
		plot_3d_pca(PC_expression, annotation, sett)
		pal = sns.cubehelix_palette(palette_size+1, start=2, rot=0, dark=0, light=0.85)
		pal = [(int(i[0]*256),int(i[1]*256),int(i[2]*256)) for i in pal]
		color_indices = map(int,pseudotime*palette_size)
		annotation["color"] = [RGBToHTMLColor(pal[i]) for i in color_indices]
		#HTML_pal = ['#a50026','#d73027','#f46d43','#fdae61','#fee08b','#ffffbf','#d9ef8b','#a6d96a','#66bd63','#1a9850','#006837']
		#annotation["color"] = [HTML_pal[i] for i in color_indices]
		
		plot_3d_pca(PC_expression, annotation, sett)
	
		#print pseudotime
	
	#plot_gene_with_pseudotime(expression_table, pseudotime.copy(), "ENST00000611179", annotation)
	
	#for tr in expression_table.columns:
	#	plot_gene_with_pseudotime(expression_table, pseudotime.copy(), tr, annotation, filename="gene_pt_plots_733/"+tr+".png")
	
	#plot_gene_with_pseudotime(expression_table, pseudotime.copy(), "ENST00000611179")

#~ cluster_dir = "/home/skevin/single_cell_pipeline/scde_input/diffex_by_trs_clusters_1_4"
def read_in_diffex(cluster_dir):

	def absoluteFilePaths(directory):
		for dirpath,_,filenames in os.walk(directory):
			for f in filenames:
				yield os.path.abspath(os.path.join(dirpath, f))

	csv_files = [x for x in absoluteFilePaths(cluster_dir) if x.endswith("_expression_values.csv")] 

	csvs = [pd.read_csv(x, delim_whitespace=True) for x in csv_files]

	csvs0 = map(lambda df: df[df.p_val < 0.05], csvs)

	itx_csvs = pd.concat(csvs0, axis=1, join="inner")
	union_csvs = pd.concat(csvs0, axis=0)

	#union_csvs.to_csv("diffex_genes_union.csv")

	csv_names = map(os.path.basename, csv_files)

	csv_names = [i.replace("_stringtie_expression_values", "") for i in csv_names]

	diffex_csvs = dict(zip(csv_names, csvs0))

	#~ for key,value in diffex_csvs.iteritems():
		#~ value.to_csv(key) 
	return diffex_csvs





if __name__ == "__main__":
	main()
