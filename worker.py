import numpy as np
import copy
from scipy.sparse import csr_matrix
import scanpy as sc
from tqdm import tqdm
import networkx as nx
from math import sqrt
from numpy import random
from numpy import linalg
import pandas as pd
import SEACells as SEACells
import tkinter

from celery import Celery
app = Celery('tasks')
app.config_from_object('celery_config')

def sparsify_assignments(A, thresh: float):
    """
    Zero out all values below a threshold in an assignment matrix
    :param A: (csr_matrix) of shape n_cells x n_SEACells containing assignment weights
    :param thresh: (float) threshold below which to zero out assignment weights
    :return: (np.array) of shape n_cells x n_SEACells containing assignment weights
    """
    A = copy.deepcopy(A)
    A[A < thresh] = 0

    # Renormalize
    A = A / A.sum(1, keepdims=True)
    A.sum(1)

    return A

# right now sparsify threshold is 0
def summarise_by_SEACell(ad, A, summarize_layer: str = 'raw'):
    """
    Construct an anndata object containing SEACell data from single cell data and assignment weights.
    :param ad: (AnnData) anndata object containing single cell data
    :param A: (csr_matrix) of shape n_cells x n_SEACells containing assignment weights
    :param summarize_layer: (str) layer of anndata
    :return: (AnnData) anndata object containing SEACell data
    """
    if summarize_layer == 'raw' and ad.raw != None:
        data = ad.raw.X
    else:
        data = ad.layers[summarize_layer]

    A = sparsify_assignments(A.T, thresh=0)
    n_cells, n_SEACells = A.shape

    print('Constructing SEACell anndata from single cells and assignment weights...')
    seacell_expressions = []
    for ix in tqdm(range(n_SEACells)):
        cell_weights = A[:, ix]
        # Construct the SEACell expression using the
        seacell_exp = data.multiply(cell_weights[:, np.newaxis]).toarray().sum(0) / cell_weights.sum()
        seacell_expressions.append(seacell_exp)

    seacell_expressions = csr_matrix(np.array(seacell_expressions))
    seacell_ad = sc.AnnData(seacell_expressions, dtype=seacell_expressions.dtype)
    seacell_ad.var_names = ad.var_names
    seacell_ad.obs['Pseudo-sizes'] = A.sum(0)
    seacell_ad.var_names = ad.var_names
    seacell_ad.obs_names = ['SEACell-' + str(i) for i in range(n_SEACells)]

    print('SEACell anndata constructed.')
    return seacell_ad

@app.task
def seacells_algorithm(fname):
    # Copy the counts to ".raw" attribute of the anndata since it is necessary for downstream analysis
    # This step should be performed after filtering

    ad = sc.read(fname)
    raw_ad = sc.AnnData(ad.X)
    raw_ad.obs_names, raw_ad.var_names = ad.obs_names, ad.var_names
    ad.raw = raw_ad

    # Normalize cells, log transform and compute highly variable genes
    sc.pp.normalize_per_cell(ad)
    sc.pp.log1p(ad)
    sc.pp.highly_variable_genes(ad, n_top_genes=1500)

    # Compute principal components -
    # Here we use 50 components. This number may also be selected by examining variance explained
    sc.tl.pca(ad, n_comps=50, use_highly_variable=True)

    ## User defined parameters

    ## Core parameters
    n_SEACells = 90
    build_kernel_on = 'X_pca' # key in ad.obsm to use for computing metacells
                            # This would be replaced by 'X_svd' for ATAC data

    ## Additional parameters
    n_waypoint_eigs = 10 # Number of eigenvalues to consider when initializing metacells

    model = SEACells.core.SEACells(ad,
                    build_kernel_on=build_kernel_on,
                    n_SEACells=n_SEACells,
                    n_waypoint_eigs=n_waypoint_eigs,
                    convergence_epsilon = 1e-5)

    
    model.construct_kernel_matrix()
    # M = model.kernel_matrix

    print('Initializing SEACells...')
    # Initialize archetypes
    model.initialize_archetypes()

    print('Running SEACells...')
    model.fit(min_iter=10, max_iter=120)

    # Compute SEACell anndata
    # SEACell_ad = summarise_by_SEACell(ad, model.A_)
    SEACell_ad = SEACells.core.summarize_by_SEACell(ad, SEACells_label='SEACell', summarize_layer='raw')

    return ad, SEACell_ad, model.A_