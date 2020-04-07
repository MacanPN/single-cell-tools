#!/usr/bin/env python
# coding: utf-8

# # import packages for velocyto

# # import packages for single-cell-tools

# In[300]:


import sys
sys.path.insert(1, '/dataVolume/storage/python_packages/single-cell-tools')
from sc_pseudotime import *


# In[302]:


get_ipython().system('mkdir data')


# # files for single-cell-tools

# In[303]:


cellset_file = "/dataVolume/storage/python_packages/single-cell-tools/resources/input/FACS_0407_2017_SHL_input_files/cell_sets_0407_SHL_20180523.csv"
settings_file = "/dataVolume/storage/python_packages/single-cell-tools/resources/input/FACS_0407_2017_SHL_input_files/plot_setting_0407_SHL_20180212.csv"

sett = settings(settings_file, cellset_file)
expression_file = "/dataVolume/storage/python_packages/single-cell-tools/resources/input/FACS_0407_2017_SHL_input_files/sunhye_census_matrix_20170407-comma_delimited.csv"
expression_table, annotation = read_expression(expression_file, sett, min_expression = 0.1, min_cells = 5)


# # reformat cell names 

# In[304]:


annotation.index = annotation.index.str.replace('X', 'shl20170407_S')


# # load loom file

# In[305]:


loom_file = "/dataVolume/storage/single_cell_projects/sc_RB_devel/20170407-SHL-FACS-Hs_proj/output/velocyto/MyTissue.loom"


# In[299]:


import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import loompy
import velocyto as vcy
import scanpy as scn
import logging
from sklearn.svm import SVR
from sklearn.linear_model import LinearRegression
from statsmodels.nonparametric.smoothers_lowess import lowess
from scipy.interpolate import interp1d

logging.basicConfig(stream=sys.stdout, format='%(asctime)s - %(levelname)s - %(message)s', level=logging.DEBUG)
get_ipython().run_line_magic('matplotlib', 'inline')
plt.rcParams['pdf.fonttype'] = 42


# # functions for velocyto

# In[301]:


# plotting utility functions
def despline():
    ax1 = plt.gca()
    # Hide the right and top spines
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    # Only show ticks on the left and bottom spines
    ax1.yaxis.set_ticks_position('left')
    ax1.xaxis.set_ticks_position('bottom')
    
def minimal_xticks(start, end):
    end_ = np.around(end, -int(np.log10(end))+1)
    xlims = np.linspace(start, end_, 5)
    xlims_tx = [""]*len(xlims)
    xlims_tx[0], xlims_tx[-1] = f"{xlims[0]:.0f}", f"{xlims[-1]:.02f}"
    plt.xticks(xlims, xlims_tx)

def minimal_yticks(start, end):
    end_ = np.around(end, -int(np.log10(end))+1)
    ylims = np.linspace(start, end_, 5)
    ylims_tx = [""]*len(ylims)
    ylims_tx[0], ylims_tx[-1] = f"{ylims[0]:.0f}", f"{ylims[-1]:.02f}"
    plt.yticks(ylims, ylims_tx)


# # Load data

# In[306]:


# Create an analysis object
vlm = vcy.VelocytoLoom(loom_file)


# In[307]:


annotation[annotation.index.isin(vlm.ca['CellID'])]


# In[308]:


annotation_index = np.isin(vlm.ca['CellID'], annotation.index)
annotation
vlm.filter_cells(annotation_index)
# annotation.index.to_numpy
# vlm.ca['CellID'].where(annotation.index)


# In[309]:


colors_dict = dict(zip(annotation.color, annotation.color))


# In[310]:


vlm.ca['CellID'].size
color_dict = {
    'default_cluster': [0.65,0.1,0.4]
}

vlm.ca['ClusterName'] = np.repeat(np.array(['default_cluster']), [vlm.ca['CellID'].size])
vlm.ca['ClusterName']
vlm.set_clusters(vlm.ca["ClusterName"], color_dict)


# # Velocity Analysis

# # filter cells by qc

# In[311]:


vlm.filter_cells(bool_array=vlm.initial_Ucell_size > np.percentile(vlm.initial_Ucell_size, 0.4))


# In[312]:


# vlm.ts = np.column_stack([vlm.ca["TSNE1"], vlm.ca["TSNE2"]])


# In[313]:


vlm.score_detection_levels(min_expr_counts=40, min_cells_express=30)
vlm.filter_genes(by_detection_levels=True)


# In[314]:


vlm.score_cv_vs_mean(200, plot=True, max_expr_avg=35)
vlm.filter_genes(by_cv_vs_mean=True)


# In[315]:


vlm.score_detection_levels(min_expr_counts=0, min_cells_express=0, min_expr_counts_U=25, min_cells_express_U=20)
vlm.score_cluster_expression(min_avg_U=0.01, min_avg_S=0.08)
vlm.filter_genes(by_detection_levels=True, by_cluster_expression=True)


# In[316]:


# best with sample and expression scaling
vlm._normalize_S(relative_size=vlm.initial_cell_size,
                 target_size=np.mean(vlm.initial_cell_size))
vlm._normalize_U(relative_size=vlm.initial_Ucell_size,
                 target_size=np.mean(vlm.initial_Ucell_size))


# In[317]:


vlm.perform_PCA()
plt.plot(np.cumsum(vlm.pca.explained_variance_ratio_)[:100])
n_comps = np.where(np.diff(np.diff(np.cumsum(vlm.pca.explained_variance_ratio_))>0.002))[0][0]
plt.axvline(n_comps, c="k")
n_comps


# In[319]:


k = 25
vlm.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True, b_sight=k*8, b_maxl=k*4, n_jobs=16)


# In[320]:


vlm.fit_gammas(limit_gamma=False, fit_offset=False)


# In[321]:


vlm.predict_U()
vlm.calculate_velocity()
vlm.calculate_shift(assumption="constant_velocity")
vlm.extrapolate_cell_at_t(delta_t=1.)


# # run tSNE

# In[322]:


from sklearn.manifold import TSNE
bh_tsne = TSNE()
vlm.ts = bh_tsne.fit_transform(vlm.pcs[:, :25])


# In[323]:


vlm.pcn = getattr(vlm, "pcs")[:,0:2]


# In[324]:


plt.figure(None, (17,2.8), dpi=80)
gs = plt.GridSpec(1,6)
for i, gn in enumerate(["ARR3","CRX"]):
    ax = plt.subplot(gs[i*3])
    try:
        ix=np.where(vlm.ra["Gene"] == gn)[0][0]
    except:
        continue
    vcy.scatter_viz(vlm.Sx_sz[ix,:], vlm.Ux_sz[ix,:], c=vlm.colorandum, s=5, alpha=0.4, rasterized=True)
    plt.title(gn)
    xnew = np.linspace(0,vlm.Sx[ix,:].max())
    plt.plot(xnew, vlm.gammas[ix] * xnew + vlm.q[ix], c="k")
    plt.ylim(0, np.max(vlm.Ux_sz[ix,:])*1.02)
    plt.xlim(0, np.max(vlm.Sx_sz[ix,:])*1.02)
    minimal_yticks(0, np.max(vlm.Ux_sz[ix,:])*1.02)
    minimal_xticks(0, np.max(vlm.Sx_sz[ix,:])*1.02)
    despline()
    
    vlm.plot_velocity_as_color(gene_name=gn, gs=gs[i*3+1], s=3, rasterized=True)

    vlm.plot_expression_as_color(gene_name=gn, gs=gs[i*3+2], s=3, rasterized=True)
    
plt.tight_layout()
plt.savefig("Fig3_selection.pdf")


# In[325]:


plt.figure(None, (16.5,15), dpi=80)
gs = plt.GridSpec(6,6)
for i, gn in enumerate(["VSX2", "OPN1SW"]):
    ax = plt.subplot(gs[i*3])
    try:
        ix=np.where(vlm.ra["Gene"] == gn)[0][0]
    except:
        continue
    vcy.scatter_viz(vlm.Sx_sz[ix,:], vlm.Ux_sz[ix,:], c=vlm.colorandum, s=5, alpha=0.4, rasterized=True)
    plt.title(gn)
    xnew = np.linspace(0,vlm.Sx[ix,:].max())
    plt.plot(xnew, vlm.gammas[ix] * xnew + vlm.q[ix], c="k")
    plt.ylim(0, np.max(vlm.Ux_sz[ix,:])*1.02)
    plt.xlim(0, np.max(vlm.Sx_sz[ix,:])*1.02)
    minimal_yticks(0, np.max(vlm.Ux_sz[ix,:])*1.02)
    minimal_xticks(0, np.max(vlm.Sx_sz[ix,:])*1.02)
    despline()
    
    vlm.plot_velocity_as_color(gene_name=gn, gs=gs[i*3+1], s=3, rasterized=True)

    vlm.plot_expression_as_color(gene_name=gn, gs=gs[i*3+2], s=3, rasterized=True)
    
plt.tight_layout()
plt.savefig("Suppl_phase_selection.pdf")


# In[326]:


vlm.estimate_transition_prob(hidim="Sx_sz", embed="pcn", transform="sqrt", psc=1,
                             n_neighbors=20, knn_random=True, sampled_fraction=0.5)


# In[328]:


vlm.calculate_embedding_shift(sigma_corr = 0.05, expression_scaling=False)
vlm.calculate_grid_arrows(smooth=0.8, steps=(40, 40), n_neighbors=30)


# In[331]:


plt.figure(None,(14,14))
quiver_scale = 60

plt.scatter(vlm.embedding[:, 0], vlm.embedding[:, 1],
            c="0.8", alpha=0.2, s=10, edgecolor="")

ix_choice = np.random.choice(vlm.embedding.shape[0], size=int(vlm.embedding.shape[0]/1.), replace=False)
plt.scatter(vlm.embedding[ix_choice, 0], vlm.embedding[ix_choice, 1],
            c="0.8", alpha=0.4, s=10, edgecolor=(0,0,0,1), lw=0.3, rasterized=True)

quiver_kwargs=dict(headaxislength=7, headlength=11, headwidth=8,linewidths=0.25, width=0.00045,edgecolors="k", color=vlm.colorandum[ix_choice], alpha=1)
plt.quiver(vlm.embedding[ix_choice, 0], vlm.embedding[ix_choice, 1],
           vlm.delta_embedding[ix_choice, 0], vlm.delta_embedding[ix_choice, 1],
           scale=quiver_scale, **quiver_kwargs)

# plt.axis("off")
plt.savefig("full_arrows.pdf")


# In[330]:


# initial divide by mean
plt.figure(None,(20,10))
vlm.plot_grid_arrows(quiver_scale=0.48,
                     scatter_kwargs_dict={"alpha":0.35, "lw":0.35, "edgecolor":"0.4", "s":38, "rasterized":True}, min_mass=24, angles='xy', scale_units='xy',
                     headaxislength=2.75, headlength=5, headwidth=4.8, minlength=1.5,
                     plot_random=True, scale_type="absolute")
plt.savefig("vectorfield.pdf")


# In[ ]:




