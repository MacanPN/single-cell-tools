#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import loompy
import velocyto as vcy
import logging
from sklearn.neighbors.kde import KernelDensity
logging.basicConfig(stream=sys.stdout, format='%(asctime)s - %(levelname)s - %(message)s', level=logging.DEBUG)
get_ipython().run_line_magic('matplotlib', 'inline')


# In[2]:


plt.rcParams['pdf.fonttype'] = 42


# In[3]:


import plotnine
import pickle 
with open ("for_seurat.pkl", "rb") as f:
    script_pkl = pickle.load(f)


# In[4]:


def rename_shl(mylist):
    mylist = [i.replace("X", "") for i in mylist]
    max_nchar = len(max(mylist, key = len))
    mylist = [i.zfill(max_nchar) for i in mylist]
    mylist = ["shl20170407-"+i for i in mylist]
    return mylist


# In[5]:


data = script_pkl["expression_table"]
pca = script_pkl["PC_expression"]
meta_data = script_pkl["annotation"]


# In[6]:


data.index = rename_shl(data.index)
pca.index = rename_shl(pca.index)
meta_data.index = rename_shl(meta_data.index)


# In[7]:


def ixs_thatsort_a2b(a: np.ndarray, b: np.ndarray, check_content: bool=True) -> np.ndarray:
    "This is super duper magic sauce to make the order of one list to be like another"
    if check_content:
        assert len(np.intersect1d(a, b)) == len(a), f"The two arrays are not matching"
    return np.argsort(a)[np.argsort(np.argsort(b))]


# In[8]:


import scvelo as scv
import anndata


# In[9]:


adata = anndata.AnnData(data, meta_data)


# In[24]:


# Crate an analysis object
adata_loom = scv.read("../20170407-SHL-FACS-Hs_proj.loom", cache = True)


# In[25]:


retained_cells = list(set(adata_loom.obs.index).intersection(set(adata.obs.index)))
retained_cells.sort()
adata_loom = adata_loom[retained_cells,:]


# In[26]:


adata_loom.var_names_make_unique()


# In[27]:


# scv.pl.proportions(adata_loom)


# #preprocess

# In[28]:


scv.pp.filter_and_normalize(adata_loom, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata_loom, n_pcs=30, n_neighbors=30)


# In[29]:


pca.shape


# In[32]:


adata_loom.obsm['X_pca'].shape
adata_loom.obsm['X_pca'][:,0:20] = pca


# In[33]:


scv.tl.velocity(adata_loom)


# In[35]:


scv.tl.velocity_graph(adata_loom)


# In[47]:


scv.pl.velocity_embedding(adata_loom, basis='pca', components='3,4')


# In[ ]:




