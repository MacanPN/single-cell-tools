#!/usr/bin/python

## @program sc-analysis
# file sc-analysis.py
# this package is built to enable temporal analysis of single cell rna-seq data

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

## what modes can be script run in
run_modes = ["find_day_correlated_pcs", "2d-pca-multiplot", "2d-pca-single", "3d-pca", "hierarchy", "pseudotime"]

## sets with parameter look like:
# operation	set_name	parameter
# for ex.: color	day_4	clue
accepted_sets_with_parameter = ["color", "outline-color", "size", "name", "shape", "superimpose", "superimpose-for-spearman"]

## sets without parameter look like:
# operation	set_name
# for ex.: remove	low_read_count_cells
accepted_sets_without_parameter = ["remove"]

## parameters supposed to be set once
accepted_parameters = ["number_of_genes"]

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
		with open(cellset_file) as f:
			for line in f:
				x = line.rstrip().split("\t")
				cell_sets[x[0]] = x[1:]
		return cell_sets
	
	def __init__(self, settings_file, cellset_file):
		self.cell_sets = self.read_cell_sets(cellset_file)
		self.sets = {}
		# inniciate all sets to empty
		for i in accepted_sets_with_parameter+accepted_sets_without_parameter:
			self.sets[i] = set()
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
			if self.run_mode not in run_modes:
				print "Unkown run mode (line 2 in settings file): ",self.run_mode
				exit(1)
			# if we're plotting pca, we want list of PCs to use
			if mode_line in ["2d-pca-single", "3d-pca"]:
				self.pcs = map(int,type_line.split("\t")[1].split(","))
				if not(
					((mode_line == "2d-pca-single")and(len(self.pcs)==2))
					or((mode_line == "3d-pca")and(len(self.pcs)==3))
					):
					print "Invalid number of PCs given! ",self.mode_line
					exit(1)
			# from third line onwards, the script reads different operations carried out on defined cell sets
			for line in f:
				if(line.startswith("#")):
					continue
				x = line.rstrip().split("\t")
				if(x[0] in accepted_sets_without_parameter):
					self.sets[x[0]].add(x[1])
				elif(x[0] in accepted_sets_with_parameter):
					self.sets[x[0]].add(tuple(x[1:3]))
				elif(x[0] in accepted_parameters):
					self.parameters[x[0]] = x[1]
				else:
					print "Unknown option:",line
					exit(1)


## function takes expression file and settings object and returns:
# - pd.DataFrame with [log transformed] expression values [genes expressed over min_expression in at least min_cells]
# - pd.DataFrame with annotations for each cell. Expression table and annotation table have the same rows
def read_expression(expression_file, settings, min_expression = 0.1, min_cells = 10, log_transform = True):
	# read expression
	expression_table = pd.read_csv(expression_file, sep=",").transpose()
	print "Read expression table with shape:",expression_table.shape
	
	# remove unwanted cells
	for s in settings.sets["remove"]:
		print "Removed cells from set:",s,settings.cell_sets[s]
		expression_table.drop(settings.cell_sets[s], inplace=True, errors="ignore")
	
	# log transform
	if(log_transform):
		expression_table += 1
		expression_table = expression_table.apply(np.log2)
		print "Log transformed data"
	# remove genes with less then min_cells expressing them
	expressed_genes = (expression_table > min_expression).sum() > min_cells
	expression_table = expression_table.loc[ : , expressed_genes]
	print "Removed genes that are not expressed >",min_expression," in at least",min_cells ,"cells"
	print "Expression table has now shape:",expression_table.shape
	
	# create annotation table and populate it with default values
	annotation = pd.DataFrame(index=expression_table.index)
	annotation["color"] = "black"
	annotation["superimpose"] = False
	annotation["superimpose-for-spearman"] = False
	annotation["size"] = 5.0
	annotation["name"] = ""
	
	for s in accepted_sets_with_parameter: # iterating over dictionary operation->set
		for i in settings.sets[s]: # iterating over set
			subset = set(settings.cell_sets[i[0]]).intersection(annotation.index)
			annotation.loc[subset,s] = i[1]
	
	annotation["size"] = pd.to_numeric(annotation["size"])
	# where outline color was not defined, set it to the color of the cell
	annotation.loc[annotation["outline-color"]!=annotation["outline-color"], "outline-color"] = annotation["color"]
	
	# define day and treatment columns
	day_labels = ["day_4","day_6","day_8","day_12"]
	treatment_labels = ["shCtrl","sh733", "sh737"]
	for i in day_labels:
		subset = set(settings.cell_sets[i]).intersection(annotation.index)
		annotation.loc[subset,"day"]=int(i.split("_")[1])
	for i in treatment_labels:
		subset = set(settings.cell_sets[i]).intersection(annotation.index)
		annotation.loc[subset,"treatment"]=i
	
	# crop annotation dataframe to only rows, that are in expression table
	annotation = annotation.loc[expression_table.index]
	return (expression_table, annotation)

## runs PCA and returns:
# - PCA transformed coordinates
# - sklearn.decomposition.pca object
def run_PCA(expression_table, annotation, n_components):
	pca = decomposition.PCA(n_components=n_components, svd_solver="full")
	expression_table_for_PCA = expression_table.loc[annotation[annotation["superimpose"]==False].index]
	print "Calculating PCA on table of shape:",expression_table_for_PCA.shape
	pca.fit(expression_table_for_PCA)
	print "Explained variance: ", pca.explained_variance_
	print "Explained variance ratio: ", pca.explained_variance_ratio_
	# transform expression using PCA vectors
	transformed_expression = pd.DataFrame(pca.transform(expression_table), index=expression_table.index, columns = range(1,n_components+1))
	return transformed_expression, pca

## save genes correlated with a PC to file
def get_isoforms_correlated_with_pc(expression_table, pc, n, filename):
	pc = int(pc)
	df = pd.Series(pca.components_[pc], index=expression_table.columns)
	df = df.reindex(df.abs().sort_values(inplace=False, ascending=False).index).iloc[0:n]
	csv_filename = settings.result_filename+"_PC"+str(pc)+".csv"
	df.to_csv(csv_filename, sep="\t")

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
			#import code; code.interact(local=locals())
			transformed_expression.loc[cells_with_this_shape].plot.scatter(
				x=pc+1,
				y=pc+2,
				ax=ax[pc/6][(pc/2)%3],
				s=ann["size"].values,
				c=ann["color"].values,
				legend=True,
				alpha=0.8,
				#edgecolor="black",
				marker = m
			)
		
		explained_variance1 = "{0:.2f}".format(pca.explained_variance_ratio_[pc]*100)+"%"
		explained_variance2 = "{0:.2f}".format(pca.explained_variance_ratio_[pc+1]*100)+"%"
		ax[pc/6][(pc/2)%3].set_xlabel("PCA "+str(pc+1)+" ["+explained_variance1+" of variance]")
		ax[pc/6][(pc/2)%3].set_ylabel("PCA "+str(pc+2)+" ["+explained_variance2+" of variance]")
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
def plot_2d_pca_single_plot(transformed_expression, annotation, pca, settings):
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
			marker = m
		)
	for cell in transformed_expression.index:
		row = transformed_expression.loc[cell,[int(settings.pcs[0]),int(settings.pcs[1])]]
		df  = transformed_expression.loc[ :  ,[int(settings.pcs[0]),int(settings.pcs[1])]]
		annotate_df(row, df, 8.0, ax)
	
	#ax.set_xlim([-100,100])
	#ax.set_ylim([-100,100])
	
	plt.xlabel("PCA "+str(settings.pcs[0]))
	plt.ylabel("PCA "+str(settings.pcs[1]))
	plt.tight_layout()
	plt.subplots_adjust(right=0.94)
	plt.savefig(settings.result_filename+"PC-"+str(settings.pcs[0])+"-"+str(settings.pcs[1])+".png", dpi=200)
	plt.show()
	#plt.close()

## create 3d PCA plot using plotly library
# arguments are:
# - pd.DataFrame with PCA transformed gene expression 
# - annotation pd.DataFrame
# - settings object
def plot_3d_pca(transformed_expression, annotation, settings, height = 1080, width = 1600):
	import plotly.plotly as py
	import plotly.graph_objs as go
	layout = dict(
		width=width,
		height=height,
		autosize=False,
		#title='Test',
		scene=dict(
			xaxis=dict(
				title="PC "+str(settings.pcs[0]),
				gridcolor='rgb(0, 0, 0)',
				zerolinecolor='rgb(255, 0, 0)',
				showbackground=True,
				backgroundcolor='#bababa'
			),
			yaxis=dict(
				title="PC "+str(settings.pcs[1]),
				gridcolor='rgb(0, 0, 0)',
				zerolinecolor='rgb(255, 0, 0)',
				showbackground=True,
				backgroundcolor='#bababa'
			),
			zaxis=dict(
				title="PC "+str(settings.pcs[2]),
				gridcolor='rgb(0, 0, 0)',
				zerolinecolor='rgb(255, 0, 0)',
				showbackground=True,
				backgroundcolor='#bababa'
			),
			aspectmode = 'manual'        
		),
	)
	data = []
	trace = dict(
		text = transformed_expression.index,# + " "+ transformed_expression["day"], #+ "\n" + transformed_expression["branch"],
		x = transformed_expression[settings.pcs[0]],
		y = transformed_expression[settings.pcs[1]],
		z = transformed_expression[settings.pcs[2]],
		type = "scatter3d",    
		mode = 'markers',
		marker = dict(
			size=annotation["size"].values,
			color=annotation["color"].values,
			symbol=annotation["shape"].values,
			line=dict(width=1) )
		)
	data.append( trace )
	fig = dict(data=data, layout=layout)
	url = py.plot(fig, filename=settings.result_filename, validate=False)

## plot hierarchycal clustering
# arguments are:
# - pd.DataFrame with PCA transformed gene expression 
# - annotation pd.DataFrame
# - settings object
# - filename for output picture
def plot_hierarchycal_clusterings(transformed_expression, annotation, settings):
	link_color = {}
	def link_color_func(node):
		return link_color[node]
	
	def colorize_links(linkage):
		l_color = {}
		n = expression_table.shape[0]
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
	
	scipy_linkage_methods = [ "complete", "average", "single", "centroid", "median", "ward"] #"single",weighted
	# plot clusterings on one magor figure
	fig,ax = plt.subplots(nrows=2, ncols=3, figsize=(50, 30))
	i=0
	for method in scipy_linkage_methods:
		linkage = sc.cluster.hierarchy.linkage(expression_table, method=method)
		link_color = colorize_links(linkage)
		dendro  = sc.cluster.hierarchy.dendrogram(
			linkage,
			ax=ax[i/3,i%3],
			labels = expression_table.index,
			link_color_func = link_color_func,
			#color_threshold = 0,
			#above_threshold_color = "black",
			count_sort = "ascending") #, title=method
		ax[i/3,i%3].set_title(method)
		tick_labels = ax[i/3,i%3].get_xmajorticklabels()
		for lbl in tick_labels:
			lbl.set_color(annotation.loc[lbl.get_text()]["color"])
		i += 1
	
	plt.tight_layout()
	plt.savefig(settings.result_filename+"-clustering.png", dpi=200)
	plt.show()

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
def find_pseudotime(transformed_expression, annotation, pca, settings):
	n_pca = len(transformed_expression.columns)
	transformed_expression["day"] = annotation["day"]
	transformed_expression_without_superimposed = transformed_expression.loc[annotation[annotation["superimpose-for-spearman"]==False].index]
	print "Finding best rotation for Spearman correlation. Shape of used table:",transformed_expression_without_superimposed.shape
	spearman = transformed_expression_without_superimposed.corr(method="spearman").loc["day",range(1,n_pca+1)].abs().sort_values(ascending=False)
	#plot_spearman correlations and explained variation
	searman_filename = settings.result_filename.replace(".png", "_spearman.png")
	width=0.4
	fig,ax = plt.subplots(figsize=(8,5))
	ax2= ax.twinx()
	spearman.plot.bar(ax=ax, width=width, position=1, color="blue")
	pd.Series(pca.explained_variance_ratio_, index=range(1,n_pca+1)).loc[spearman.index].plot.bar(ax=ax2, width=width, position=0, color="red")
	ax.set_xlabel("PC component")
	ax.set_ylabel("Spearman correlation\nto days [blue]")
	ax2.set_ylabel("% variance explained [red]")
	plt.tight_layout()
	low,high = plt.xlim()
	plt.xlim(low-0.5, high)
	plt.savefig(searman_filename, dpi=200)
	settings.pcs = map(str,spearman.iloc[0:2].index)
	
	# find best rotation
	best_angle = 0
	best_spearman = 0
	for a in range(0,360):
		te = rotate_expression(transformed_expression_without_superimposed, int(settings.pcs[0]), int(settings.pcs[1]), a)
		spearman = te.corr(method="spearman").loc["day",int(settings.pcs[0])]
		#print "Trying angle: ",a," spearman: ",spearman
		if(spearman > best_spearman):
			best_angle = a
			best_spearman = spearman
		
	del(transformed_expression["day"])
	print settings.pcs
	print "Best rotation: ",best_angle
	rotated_expression = rotate_expression(transformed_expression, int(settings.pcs[0]), int(settings.pcs[1]), best_angle)
	#plot_2d_pca_single_plot(rotated_expression, filename = settings.result_filename.replace(".png","-rotated.png"))
	pt = rotated_expression[int(settings.pcs[0])]
	pt.name = "pseudotime"
	return pt

## plots gene expression over pseudotime
# arguments are:
# - pd.DataFrame with gene expression 
# - pd.Series with pseudotime coordinates for each cell
# - Ensamble transcript ID
def plot_gene_with_pseudotime(exp, pseudotime, transcript_id):
	expression_over_pseudotime = pd.DataFrame(pseudotime)
	expression_over_pseudotime["expression"] = exp[transcript_id]
	expression_over_pseudotime.plot.scatter(x="pseudotime", y="expression")
	plt.show()
	

## main function
#  when run separately, program expects following arguments:
# - argv[1] = comma separated file with expression values
# - argv[2] = file with cell sets (see settings.read_cell_sets())
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
		plot_hierarchycal_clusterings(PC_expression, annotation, sett)
	elif(sett.run_mode=="pseudotime"):
		pseudotime = find_pseudotime(PC_expression, annotation, pca, sett)
	
	#plot_gene_with_pseudotime(expression_table, pseudotime.copy(), "ENST00000611179")

if __name__ == "__main__":
	main()
