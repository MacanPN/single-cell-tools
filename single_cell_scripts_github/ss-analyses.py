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

# get parameters
expression_file = sys.argv[1]
cellset_file    = sys.argv[2]
settings_file   = sys.argv[3]
n_pca = 20

# cell sets
cell_sets = {}
with open(cellset_file) as f:
	for line in f:
		x = line.rstrip().split("\t")
		cell_sets[x[0]] = x[1:]

class settings:
	def __init__(self, settings_file):
		self.remove_sets = []
		self.set_colors  = []
		self.set_outline_colors  = []
		self.set_sizes   = []
		self.set_names   = []
		self.set_shapes  = []
		self.set_superimpose = []
		self.set_superimpose_spearman = []
		self.pcs = []
		self.plot_type = ""
		self.number_of_genes = 1000
		with open(settings_file) as f:
			self.result_filename = "PCA_plots/"+f.readline().rstrip()+".png"
			valid_plot_types = ["find_day_correlated_pcs", "2d-pca-multiplot", "2d-pca-single", "3d-pca", "hierarchy"]
			type_line = f.readline().rstrip()
			self.plot_type = type_line.split("\t")[0]
			if self.plot_type in ["2d-pca-single", "3d-pca"]:
				self.pcs = type_line.split("\t")[1].split(",")
			for line in f:
				x = line.rstrip().split("\t")
				if(x[0] == "remove"):
					self.remove_sets.append(x[1])
				elif(x[0] == "color"):
					self.set_colors.append(x[1:3])
				elif(x[0] == "outline-color"):
					self.set_outline_colors.append(x[1:3])
				elif(x[0] == "size"):
					self.set_sizes.append(x[1:3])
				elif(x[0] == "name"):
					self.set_names.append(x[1:3])
				elif(x[0] == "shape"):
					self.set_shapes.append(x[1:3])
				elif(x[0] == "superimpose"):
					self.set_superimpose.append(x[1:3])
				elif(x[0] == "superimpose-for-spearman"):
					self.set_superimpose_spearman.append(x[1:3])
				elif(x[0] == "number_of_genes"):
					self.number_of_genes = int(x[1])
				else:
					print "Unknown option:",line
				#exit(1)
settings = settings(settings_file)
print settings.plot_type
# read expression
expression_table = pd.read_csv(expression_file, sep=",").transpose()
print "Read expression table with shape:",expression_table.shape
# create annotation table
annotation = pd.DataFrame(index=expression_table.index)
annotation["color"] = "black"
annotation["superimpose"] = False
annotation["superimpose-for-spearman"] = False
annotation["size"] = 5.0
annotation["name"] = ""

for i in settings.set_colors:
	annotation.loc[cell_sets[i[0]],"color"] = i[1]

annotation["outline-color"] = annotation["color"].copy()
for i in settings.set_outline_colors:
	annotation.loc[cell_sets[i[0]],"outline-color"] = i[1]
for i in settings.set_sizes:
	annotation.loc[cell_sets[i[0]],"size"] = float(i[1])
for i in settings.set_names:
	annotation.loc[cell_sets[i[0]],"name"] = i[1]
for i in settings.set_shapes:
	annotation.loc[cell_sets[i[0]],"shape"] = i[1]
for i in settings.set_superimpose:
	annotation.loc[cell_sets[i[0]],"superimpose"] = True
for i in settings.set_superimpose_spearman:
	annotation.loc[cell_sets[i[0]],"superimpose-for-spearman"] = True
day_labels = ["day_4","day_6","day_8","day_12"]
treatment_labels = ["shCtrl","sh733", "sh737"]
for i in day_labels:
	annotation.loc[cell_sets[i],"day"]=int(i.split("_")[1])
for i in treatment_labels:
	annotation.loc[cell_sets[i],"treatment"]=i

# log transform
expression_table += 1
expression_table = expression_table.apply(np.log2)
print "Log transformed data"
# remove genes with <10 cells expressing them
min_expression = 0.1
min_cells = 10
expressed_genes = (expression_table > min_expression).sum() > min_cells
expression_table = expression_table.loc[ : , expressed_genes]
print "Removed genes that are not expressed >0.1 in at least 10 cells"
# remove unwanted cells
for s in settings.remove_sets:
	print "Removed cells from set:",s,cell_sets[s]
	expression_table.drop(cell_sets[s], inplace=True, errors="ignore")
print "expression table has now shape:",expression_table.shape

annotation = annotation.loc[expression_table.index]

# run PCA
pca = decomposition.PCA(n_components=n_pca, svd_solver="full")
expression_table_for_PCA = expression_table.loc[annotation[annotation["superimpose"]==False].index]
print "Calculating PCA on table of shape:",expression_table_for_PCA.shape
pca.fit(expression_table_for_PCA)
print "Explained variance: ", pca.explained_variance_
print "Explained variance ratio: ", pca.explained_variance_ratio_
# transform expression using PCA vectors
transformed_expression = pd.DataFrame(pca.transform(expression_table), index=expression_table.index, columns = range(1,n_pca+1))

# save genes correlated with PCs to file(s)
def get_isoforms_correlated_with_pc(expression_table, pc, n, filename=None):
	if(filename == None):
		filename = settings.result_filename
	pc = int(pc)
	df = pd.Series(pca.components_[pc], index=expression_table.columns)
	df = df.reindex(df.abs().sort_values(inplace=False, ascending=False).index).iloc[0:n]
	csv_filename = settings.result_filename+"_PC"+str(pc)+".csv"
	df.to_csv(csv_filename, sep="\t")

# annotate point on axis if it's far enough
def annotate_df(row,df,min_dist,ax):
		dist = (df - row).abs().sum(axis=1).sort_values()[1]
		#print dist
		if(dist > min_dist):
			ax.annotate(row.name, list(row.values),
				xytext=(5,-3), 
				textcoords='offset points',
				size=10, 
				color='darkslategrey')

# pandas 2d plots:
def plot_2d_pca_multiplot(transformed_expression):
	fig, ax = plt.subplots(2,3, figsize=(15,10))
	markers = list(annotation["shape"].unique())
	for pc in range(0,12,2):
		for m in markers:
			cells_with_this_shape = annotation["shape"]==m
			ann = annotation.loc[cells_with_this_shape]
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
	plt.savefig(settings.result_filename, dpi=200)
	plt.show()

def plot_2d_pca_single_plot(transformed_expression, filename=None):
	if(filename == None):
		filename = settings.result_filename
	get_isoforms_correlated_with_pc(expression_table, settings.pcs[0],settings.number_of_genes)
	get_isoforms_correlated_with_pc(expression_table, settings.pcs[1],settings.number_of_genes)
	fig,ax = plt.subplots(figsize=(5,5))
	markers = list(annotation["shape"].unique())
	for m in markers:
		cells_with_this_shape = annotation["shape"]==m
		ann = annotation.loc[cells_with_this_shape]
		transformed_expression.loc[cells_with_this_shape].plot.scatter(
			x=int(settings.pcs[0]),
			y=int(settings.pcs[1]),
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
	
	plt.xlabel("PCA "+settings.pcs[0])
	plt.ylabel("PCA "+settings.pcs[1])
	plt.tight_layout()
	plt.subplots_adjust(right=0.94)
	plt.savefig(filename, dpi=200)
	#plt.show()
	plt.close()

# plotly 3d plot
def plot_using_plotly(transformed_expression,pc):
	get_isoforms_correlated_with_pc(expression_table, settings.pcs[0],settings.number_of_genes)
	get_isoforms_correlated_with_pc(expression_table, settings.pcs[1],settings.number_of_genes)
	get_isoforms_correlated_with_pc(expression_table, settings.pcs[2],settings.number_of_genes)
	import plotly.plotly as py
	import plotly.graph_objs as go
	layout = dict(
		width=1600,
		height=1080,
		autosize=False,
		#title='Test',
		scene=dict(
			xaxis=dict(
				title="PC "+str(pc[0]),
				gridcolor='rgb(0, 0, 0)',
				zerolinecolor='rgb(255, 0, 0)',
				showbackground=True,
				backgroundcolor='#bababa'
			),
			yaxis=dict(
				title="PC "+str(pc[1]),
				gridcolor='rgb(0, 0, 0)',
				zerolinecolor='rgb(255, 0, 0)',
				showbackground=True,
				backgroundcolor='#bababa'
			),
			zaxis=dict(
				title="PC "+str(pc[2]),
				gridcolor='rgb(0, 0, 0)',
				zerolinecolor='rgb(255, 0, 0)',
				showbackground=True,
				backgroundcolor='#bababa'
			),
			#aspectratio = dict( x=1, y=1, z=0.7 ),
			aspectmode = 'manual'        
		),
	)
	data = []
	trace = dict(
		text = transformed_expression.index,# + " "+ transformed_expression["day"], #+ "\n" + transformed_expression["branch"],
		x = transformed_expression[pc[0]],
		y = transformed_expression[pc[1]],
		z = transformed_expression[pc[2]],
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
	url = py.plot(fig, filename='single_cell-3d-pca', validate=False)

# plot hierarchycal clustering
def plot_hierarchycal_clusterings(transformed_expression):
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
	plt.savefig(settings.result_filename, dpi=200)

def rotate_expression(transformed_expression,x,y,angle):
	theta = math.radians(angle)
	ret = transformed_expression.copy()
	ret[x] = transformed_expression[x]*math.cos(theta) - transformed_expression[y]*math.sin(theta)
	ret[y] = transformed_expression[x]*math.sin(theta) + transformed_expression[y]*math.cos(theta)
	return ret

# find most correlated pcs
def find_pseudotime(transformed_expression):
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
		print "Trying angle: ",a," spearman: ",spearman
		if(spearman > best_spearman):
			best_angle = a
			best_spearman = spearman
		
	del(transformed_expression["day"])
	print settings.pcs
	print "Best rotation: ",best_angle
	plot_2d_pca_single_plot(transformed_expression, filename = settings.result_filename)
	rotated_expression = rotate_expression(transformed_expression, int(settings.pcs[0]), int(settings.pcs[1]), best_angle)
	plot_2d_pca_single_plot(rotated_expression, filename = settings.result_filename.replace(".png","-rotated.png"))
	pt = rotated_expression[int(settings.pcs[0])]
	pt.name = "pseudotime"
	return pt

def plot_gene_with_pseudotime(exp, pseudotime, gene):
	expression_over_pseudotime = pd.DataFrame(pseudotime)
	expression_over_pseudotime["expression"] = exp[gene]
	expression_over_pseudotime.plot.scatter(x="pseudotime", y="expression")
	plt.show()
	
print "Now plotting:",settings.plot_type

if(settings.plot_type=="2d-pca-multiplot"):
	plot_2d_pca_multiplot(transformed_expression)
elif(settings.plot_type=="2d-pca-single"):
	plot_2d_pca_single_plot(transformed_expression)
elif(settings.plot_type=="3d-pca"):
	plot_using_plotly(transformed_expression,map(int,settings.pcs))
elif(settings.plot_type=="hierarchy"):
	plot_hierarchycal_clusterings(transformed_expression)
elif(settings.plot_type=="pseudotime"):
	pseudotime = find_pseudotime(transformed_expression)

#plot_gene_with_pseudotime(expression_table, pseudotime.copy(), "ENST00000611179")

