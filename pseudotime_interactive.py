#!/usr/bin/env python

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
import chart_studio.plotly
import colorsys
import plotly.express as px

#~ import ipdb

from sc_pseudotime import *

parser = argparse.ArgumentParser(description="runs pseudotime_interactive")
parser.add_argument("-e", "--expression-matrix", dest="expr_mat", default="resources/dshayler_input/3_Fetal_Seq/fetal_census_matrix.csv", help="gene by cell matrix of expression values", metavar="EXPR")
parser.add_argument("-c", "--cell-sets", dest="cell_sets", default="resources/dshayler_input/3_Fetal_Seq/ThreeSeq_Fetal_Metadata_101718.csv", help="cell sets", metavar="CELL_SETS")
parser.add_argument("-p", "--plot-settings", dest="plot_settings", default="resources/dshayler_input/3_Fetal_Seq/101718_3d_PCA_No_Bad_Reads_NoVSX2_3FR_Seq_grp.txt", help="plot settings", metavar="PLOT_SETTINGS")
parser.add_argument("-n", "--session-name", dest="session_name", help="a name to give to this analysis session for reproducbility", metavar="SESSION_NAME", required=False)


try:
    options = parser.parse_args()
except SystemExit as err: 
    if err.code == 2: 
        parser.print_help()
        sys.exit(0)
 
# ~ load datasets
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

def assign_clusters_using_hierarch(colnm=None, colval=None):
    pc_set = "Which PCs would you like to use for clustering? [type comma separated list, list can also include ranges 1-5,8] "
    scipy_linkage_methods = [ "complete", "average", "single", "centroid", "median", "ward"]
    cluster_on_pcs = list_from_ranges(input(pc_set))
    number_of_clusters = int(input("How many clusters would you like to generate? "))
    sett.num_clusters = number_of_clusters
    method = input("Which clustering would you like to use: "+", ".join(scipy_linkage_methods)+": ")
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
            time = float(input("Assign time for cluster shown in "+cluster_colors[i]+": "))
            clusters.append( (time,subset_annotation.loc[subset_annotation["color"]==cluster_colors[i]].index) )
        clusters.sort(key=lambda by_first: by_first[0])
        dendro = plot_hierarchical_clustering(subset_PC_expression[cluster_on_pcs], subset_annotation, method=method, sett=sett)
    else:
        linkage = sc.cluster.hierarchy.linkage(PC_expression[cluster_on_pcs], method=method)
        clusters_without_time = get_cluster_labels(linkage, number_of_clusters, PC_expression.index)
        cluster_colors = ["blue", "red", "orange", "purple", "green", "brown", "black", "gray", "lawngreen", "magenta", "lightpink", "indigo", "lightblue", "lightgoldenrod1", "mediumpurple2"]
        print("Now plotting clusters")
        subset_annotation = change_annotation_colors_to_clusters(clusters_without_time, annotation, cluster_colors)
        clusters = []
        sett.pcs = cluster_on_pcs[:3]
        plot_3d_pca(PC_expression, subset_annotation, sett)
        for i in range(0,number_of_clusters):
            time = float(input("Assign time for cluster shown in "+cluster_colors[i]+": "))
            clusters.append( (time,subset_annotation.loc[subset_annotation["color"]==cluster_colors[i]].index) )
        clusters.sort(key=lambda by_first: by_first[0])
        dendro = plot_hierarchical_clustering(PC_expression[cluster_on_pcs], subset_annotation, method=method, sett=sett)
    sett.subset = 'param'
    return clusters, dendro
    
def assign_clusters_using_file(cluster_file):
    sett.append("cluster", cluster_file)
    sett.n_clusters = len(sett.sets['cluster'])
    cluster_colors = ["blue", "red", "orange", "purple", "green", "brown", "black", "gray", "lawngreen", "magenta", "lightpink", "indigo", "lightblue", "lightgoldenrod1", "mediumpurple2"]
    clusters = {}
    number_of_clusters = sett.n_clusters
    for i in range(number_of_clusters):
        clusters[i] = set(sett.cell_sets[sett.sets['cluster'][i]])
    subset_annotation, subset_PC_expression = subset_pc_by_clusters(PC_expression, clusters)
    subset_annotation = change_annotation_colors_to_clusters(clusters, subset_annotation, cluster_colors)
    clusters = []
    plot_3d_pca(subset_PC_expression, subset_annotation, sett)
    for i in range(0,number_of_clusters):
        time = float(input("Assign time for cluster shown in "+cluster_colors[i]+": "))
        clusters.append( (time,subset_annotation.loc[subset_annotation["color"]==cluster_colors[i]].index) )
    clusters.sort(key=lambda by_first: by_first[0])
    # remove empty clusters before proceeding
    filt_clusters = []
    for i in clusters:
        if(len(i[1]) > 0):
            filt_clusters.append(i)
    # set subset flag
    sett.subset = 'cluster'
    return subset_annotation, subset_PC_expression, filt_clusters

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
        
def print_dendro(dendros, method):
    oldd = dict(zip(dendros[method]['ivl'], dendros[method]['color_list']))  

    newd = {}
    for i,d in oldd.items():
        newd.setdefault(d,[]).append(i)
    print(newd)
    return(newd)

def retrieve_subset_param():
    colnm = input("What metadata should be used to subset the data? (ex. treatment, age, etc.) ")
    colval = input("What values should be used to subset the data? (ex. shCtrl, sh842,). Providing no value will prevent subsetting ").split(",")
    if colval != "":
      sett.subset = "param"
    return colnm, colval

def subset_pc_by_param(pc_expression, colnm, colval):
    if not all(colval):
        # ~ clusters_without_time = get_cluster_labels(linkage, number_of_clusters, subset_PC_expression.index)
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
        
def subset_pc_by_clusters(pc_expression, clusters):
    cluster_cells = set.union(*clusters.values())
    subset_annotation = annotation[annotation.index.isin(cluster_cells)]   
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
    pcs = map(int,input("Which PCs would you like on the plot? (type comma separated list, such as 1,3,4) ").split(","))
    sett.pcs = pcs
    
    ctrl_annotation, ctrl_pc_expression = subset_pc_by_param(PC_expression, colnm, colval)
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
    [H]    Plot Hierarchical Clustering
    [P]    Plot PCA
    [L]    Assign clusters according to time Labels (like day_4 ... )
    [C]    Assign clusters using hierarchical clustering
    [U]    Assign clusters using file
    [D]    Find Most Correlated and Most Discriminating (treat v ctrl) PCs
    [N]    Normalize centroids
    [G]    Plot PCA Colored by Expression of Marker Genes
    [S]    Calculate pseudotime for cells using times assigned to clusters
    [O]    Output clusters (so they can be copied to a file)
    [F]    Save generated pseudotime to file
    [M]    Save features correlated with pc to file
    [T]    Run tSNE
    [X]    Exit
    """
    action = input(question).upper()
    if(action == "X"):
        break
        #~ exit()
    elif(action == "H"):
        
        color_scheme = input("How would you like to color dendrogram? ('retain' to keep current pca plot colors, 'overwrite' to use new colors) ")
        if color_scheme == "overwrite":
          number_of_clusters = int(input("How many clusters would you like to generate? "))
        else:
          number_of_clusters = 3
        sett.num_clusters = number_of_clusters
        print("plotting...\n dendrogram will be saved as a .pdf shortly")
        try:
            subset_annotation
        except:
            dendros = plot_all_hierarchical_clusterings(PC_expression, annotation, color_scheme, sett)
        else:
            dendros = plot_all_hierarchical_clusterings(PC_expression, subset_annotation, color_scheme, sett)
    elif(action == "P"):
        colnm, colvalp = retrieve_subset_param()
        pcs = [int(i) for i in input("Which PCs would you like on the plot? (type comma separated list, such as 1,3,4) ").split(",")]
        
        sett.pcs = pcs
        print("plotting...\n the plot will open in your web browser shortly")
        # IPython.embed()
        # If using settings as specified in initial settings files
        if (sett.subset == 'None'):        
        # ~ if not all(colvalp):
            fig = plot_3d_pca(PC_expression, annotation, sett, clusters = clusters)
        elif (sett.subset == 'param'):
            subset_annotation, subset_PC_expression = subset_pc_by_param(PC_expression, colnm, colvalp)
            plot_3d_pca(subset_PC_expression, subset_annotation, sett, clusters = clusters)
            del subset_annotation, subset_PC_expression
        # elif (sett.subset == 'param'):
        #   
        #     if (colvalp == colval):
        #         plot_3d_pca(subset_PC_expression, subset_annotation, sett, clusters = subset_clusters)
        #     elif set(colvalp).issubset(colval):
        #         #~ ipdb.set_trace()
        #         old_colors = subset_annotation["color"]
        #         
        #         subset_annotation, subset_PC_expression = subset_pc_by_param(subset_PC_expression, colnm, colvalp)
        #         subset_annotation.loc[:,"color"] = old_colors[old_colors.index.isin(subset_annotation.index)]
        #         new_clusters = [(i, c[c.isin(subset_annotation.index)]) for i,c in subset_clusters]
        #         plot_3d_pca(subset_PC_expression, subset_annotation, sett, clusters = new_clusters)
        #         del subset_annotation, subset_PC_expression
        #     else:
        #         subset_annotation, subset_PC_expression = subset_pc_by_param(PC_expression, colnm, colvalp)
        #         plot_3d_pca(subset_PC_expression, subset_annotation, sett, clusters = clusters)
        #         del subset_annotation, subset_PC_expression
        elif (sett.subset == 'cluster'):
            plot_3d_pca(subset_PC_expression, subset_annotation, sett, clusters = subset_clusters)
            
    elif(action == "L"):
        colnm, colval = retrieve_subset_param()
        subset_annotation, subset_PC_expression = subset_pc_by_param(PC_expression, colnm, colval)
        subset_clusters = time_clusters_from_annotations(subset_annotation)
        pcs = [int(i) for i in input("Which PCs would you like on the plot? (type comma separated list, such as 1,3,4) ").split(",")]
        print("Time clusters were assigned according to labels")
    
    elif(action == "D"):
        colnm, colval = retrieve_subset_param()
        subset_annotation, subset_PC_expression = subset_pc_by_param(PC_expression, colnm, colval)
        # ~ pcs = map(int,input("Which PCs would you like on the plot? (type comma separated list, such as 1,3,4) ").split(","))
        find_pseudotime(subset_PC_expression, subset_annotation, pca, sett)
        print("Showing PCS most correlated with time")
        
    elif(action == "C"):
        colnm, colval = retrieve_subset_param()
        subset_annotation, subset_PC_expression = subset_pc_by_param(PC_expression, colnm, colval)
        subset_clusters, dendro = assign_clusters_using_hierarch(colnm, colval)
        print("Time clusters were assigned according to hierarchical clustering")
        filename = input("Enter file name you'd like to save clustering plot as (preferably ending with .pdf) ")
        plt.savefig(filename)
 
    elif(action == "U"):
        cluster_file = input("Provide path to file with cluster info (in cell_settings format) ")
        # ~ cluster_file = "runtime_settings.csv"
        subset_annotation, subset_PC_expression, subset_clusters = assign_clusters_using_file(cluster_file)
        print("Time clusters were assigned according to specified file ")
    
    elif(action == "N"):
        test = normalize_centroids(subset_pc_expression)
        url = plotly.offline.plot(test, filename="resources/normalize_centroids.html", validate=False, auto_open=False)
        print(url)
        
    elif(action == "S"):
        colnm, colval = retrieve_subset_param()
        subset_annotation, subset_PC_expression = subset_pc_by_param(PC_expression, colnm, colval)
        pseudotime, centroids = calculate_pseudotime_using_cluster_times(subset_PC_expression, subset_annotation, subset_clusters, sett)
        
    elif(action == "O"):
        cluster_source = input("retrieve clusters from 3d or all dimension clustering (3d, all)? ")
        if cluster_source == "3d":
            print_clusters(subset_clusters)
        elif cluster_source == "all":
            cluster_method = input("which hierarch. clustering alg. (complete, average, ward, etc.)? ")
            print_dendro(dendros, cluster_method)
    
    elif(action == "G"):
        colnm, colval = retrieve_subset_param()
        subset_annotation, subset_PC_expression = subset_pc_by_param(PC_expression, colnm, colval)
        pcs = [int(i) for i in input("Which PCs would you like on the plot? (type comma separated list, such as 1,3,4) ").split(",")]
        sett.pcs = pcs
        # bins = int(input("How many bins would you like to quantile? "))
        bins = 5
        sett.bins = bins
        features = input("Which genes would you like to plot (type comma separated list, such as RB1,RXRG,ARR3) ").split(",")
        feat_type = input("Plot by gene (g) or by transcript (t)? ")
        
        #palette_size = bins
        #pal = sns.color_palette("coolwarm", palette_size+1)
        #pal = [(int(i[0]*256),int(i[1]*256),int(i[2]*256)) for i in pal]
        bin_colors = ["grey", "blue", "#21f2e4", "orange", "red", "#51040a"]
        bin_col_dict = dict(zip(range(0,bins), bin_colors))
	#bin_colors = cl.to_rgb(pal)
        #bin_col_dict = dict(zip(range(0,bins), bin_colors))
        if feat_type == "g":
            for i in features:
                plot_3d_pca(subset_PC_expression, subset_annotation, sett, clusters = clusters, features=i, bin_col_dict=bin_col_dict, expression_table=expression_table, feat_type = feat_type)
        elif feat_type == "t":
            for i in features:
                sett.result_filename = sett.result_filename+"_"+i
                markers = list(annotation["shape"].unique())
                mg = mygene.MyGeneInfo()
                # ~ for i in enumerate(features): 
                gene_info = mg.querymany(i, scopes='symbol', fields='ensembl.transcript')[0]
                if len(gene_info['ensembl']) > 1:
                    trx = gene_info['ensembl'][0]['transcript']
                else:
                    trx = gene_info['ensembl']['transcript']
                for j in trx:
                    if j in expression_table.columns:
                        plot_3d_pca(subset_PC_expression, subset_annotation, sett, clusters = clusters, features=j, bin_col_dict=bin_col_dict, expression_table=expression_table, feat_type = feat_type)
        # ~ plot_marker_gene_quantile(expression_table, subset_PC_expression, subset_annotation, pcs, sett, marker_genes)
        
        
    elif(action == "M"):
        filename = input("Enter file name you'd like to save correlated features to: ")
        top_n = input("How many genes from each pc would you like to save?; if blank will use "+sett.parameters['number_of_genes']+" (from plot_settings) ")
        if top_n == "":
            top_n = int(sett.parameters['number_of_genes'])
        else:
            top_n = int(top_n)
        pcs = [int(i) for i in input("Which PCs would you like on the plot? (type comma separated list, such as 1,3,4) ").split(",")]
        pc_corr_trs = get_isoforms_correlated_pc_set(pca, expression_table, pcs, top_n, filename)
        csv_filename = filename+"_pcs_"+"_".join(map(str, pcs))+".csv"
        print("saving as "+csv_filename)
        
    elif(action=="F"):
        if("pseudotime" not in globals()):
            print("Pseudotime was not yet generated!")
            continue
        filename = input("Enter file name you'd like to save pseudotime as (preferably ending with .csv) ")
        pseudotime.to_csv(filename, sep="\t")
    elif(action=="T"):
        
        # ~ for learning_rate in [10,50,100,200,400,1000]:
            # ~ for perplexity in [5,8,10,20,30,50]:
                # ~ print "fitting tSNE, step:",step
                # ~ tsne = sklearn.manifold.TSNE(
                    # ~ n_components = 2,
                    # ~ learning_rate = learning_rate,
                    # ~ perplexity = perplexity,
                    # ~ n_iter = 10000,
                    # ~ n_iter_without_progress = 1000
                # ~ )
                # ~ tsne_transformed_expression_2d = pd.DataFrame(tsne.fit_transform(PC_expression.values), index=PC_expression.index, columns=["x","y"])
                # ~ fig,ax = plt.subplots(figsize=(15, 10))
                # ~ ax = tsne_transformed_expression_2d.plot.scatter(x="x", y="y", s = annotation["day"].values, c=annotation["treatment"].values, ax=ax)
                # ~ ab=tsne_transformed_expression_2d.apply(annotate_df, axis=1)
                # ~ fn = "PCA_tSNE_training_lr_per/tSNE_lr."+str(learning_rate)+"_perplexity."+str(perplexity)+".png"
                # ~ plt.savefig(fn)
                # ~ plt.close()

        tsne3d = sklearn.manifold.TSNE(n_components=3, learning_rate=20, n_iter=10000, n_iter_without_progress=1000, perplexity = 10)
        tsne_transformed_expression_3d = pd.DataFrame(tsne3d.fit_transform(PC_expression.values), index=PC_expression.index, columns=["x","y","z"])
        comb = pd.concat([tsne_transformed_expression_3d, annotation], axis=1)
        # plotly 3d plot
        def plot_using_plotly(transformed_expression):
            import plotly.plotly as py
            import plotly.graph_objs as go
            layout = dict(
            width=1600,
            height=1080,
            autosize=False,
            #title='Test',
            scene=dict(
                xaxis=dict(
                    gridcolor='rgb(0, 0, 0)',
                    zerolinecolor='rgb(255, 0, 0)',
                    showbackground=True,
                    backgroundcolor='#bababa'
                ),
                yaxis=dict(
                    gridcolor='rgb(0, 0, 0)',
                    zerolinecolor='rgb(255, 0, 0)',
                    showbackground=True,
                    backgroundcolor='#bababa'
                ),
                zaxis=dict(
                    #title="PC "+str(pc[2]+1),
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
            traces = comb["name"].unique()
            for t in traces:
              
                trace = dict(
                    text = transformed_expression.index, #+ "\n" + transformed_expression["branch"],
                    x = transformed_expression["x"],
                    y = transformed_expression["y"],
                    z = transformed_expression["z"],
                    type = "scatter3d",    
                    mode = 'markers',
                    opacity = 0.80,
                    marker = dict(
                    size=comb.loc[comb["name"]==t,"size"].values,
                    # ~ color=trx_df.color,
                    color=comb.loc[comb["name"]==t,"color"].values,
                    symbol=comb.loc[comb["name"]==t,"shape"].apply(shape_matplotlib2plotly).values,
                    line=dict(width=1) )
                )
                
                data.append(trace)
            fig = dict(data=data, layout=layout)
            url = plotly.offline.plot(fig, filename='resources/single_cell-3d-tSNE', validate=False, auto_open=False)
            
            plot_using_plotly(tsne_transformed_expression_3d)
        
    elif(action=="Q"):
        cluster_dir = input("Enter location of diffex csvs ")
        diffex_csvs = read_in_diffex(cluster_dir)
        plot_heatmap(expression_table, annotation, dendro)
        
    elif(action == "I"):
        IPython.embed()
    




#HTML_pal = ['#a50026','#d73027','#f46d43','#fdae61','#fee08b','#ffffbf','#d9ef8b','#a6d96a','#66bd63','#1a9850','#006837']
#IPython.embed()
#annotation["color"] = [HTML_pal[i] for i in color_indices]
