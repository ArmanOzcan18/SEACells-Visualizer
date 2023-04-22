import numpy as np
import copy
from scipy.sparse import csr_matrix
import scanpy as sc
from tqdm import tqdm
import networkx as nx
from math import sqrt


def triangle_info(ad, A, seacell_1 = None, seacell_2 = None, seacell_3 = None):
    """
    Given three SEACells, return the adjacency matrix, labels, and coordinates for the triangle that they form with the cells inside.
    """
    v1 = [0,sqrt(3)]
    v2 = [1,0]
    v3 = [-1,0]
    tri_coords = [v1,v2,v3]
    strength = [0, 0, 0]

    if(seacell_1 == None or seacell_2 == None or seacell_3 == None):
        tri_adj_matrix = np.zeros(shape=(len(tri_coords), len(tri_coords)), dtype=np.int8)
        tri_adj_matrix[0, 1:3], tri_adj_matrix[1:3, 0], tri_adj_matrix[2, 1], tri_adj_matrix[1, 2] = 1, 1, 1, 1
        tri_labels = ["SEACell", "SEACell", "SEACell"]
        return tri_adj_matrix, tri_labels, tri_coords, strength

    tri_labels = ["SEACell-{}".format(seacell_1), "SEACell-{}".format(seacell_2), "SEACell-{}".format(seacell_3)]
    for i in range(A.shape[1]):
        weight_1 = A[seacell_1, i]
        weight_2 = A[seacell_2, i]
        weight_3 = A[seacell_3, i]
        sum = weight_1 + weight_2 + weight_3
        if(sum < 0.05):
            continue
        weight_1, weight_2, weight_3 = weight_1/sum, weight_2/sum, weight_3/sum
        coordinate = [weight_2 - weight_3, weight_1 * sqrt(3)]
        tri_coords.append(coordinate)
        strength.append(int(sum/0.05)/20)
        tri_labels.append(ad.obs.index[i])

    tri_adj_matrix = np.zeros(shape=(len(tri_coords), len(tri_coords)), dtype=np.int8)
    tri_adj_matrix[0, 1:3], tri_adj_matrix[1:3, 0], tri_adj_matrix[2, 1], tri_adj_matrix[1, 2] = 1, 1, 1, 1
    return tri_adj_matrix, tri_labels, tri_coords, strength

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
    #why is this ordered?
    seacell_ad.obs_names = ['SEACell-' + str(i) for i in range(n_SEACells)]

    print('SEACell anndata constructed.')
    return seacell_ad


def seacells_by_weights(A, threshold):
    """
    Return a list of proper SEACells assignments for each cell by thresholding the weights.
    """
    count_of_one_assignments= 0
    count_of_two_assignments = 0
    list_of_SEACell_assignments_per_cell = []
    for i in range(A.shape[1]):
        if (max(A[:, i]) == 1):
            count_of_one_assignments += 1
        no = len(A[:, i][A[:, i] > threshold])
        if (no < 3):
            count_of_two_assignments += 1
            continue
        ind = np.argsort(A[:, i])[::-1][:no]
        list_of_SEACell_assignments_per_cell.append(ind)
    return list_of_SEACell_assignments_per_cell

def compute_nn_triangles(adjacency_list):
    """
    Compute the triangles in the SEACElls graph using the adjacency list that the nearest neighbor method gave us.
    """
    nn_triangles = []
    in_triangles = set()
    data = [set() for x in range(len(adjacency_list))]
    for s in range(len(adjacency_list)):
        for t in adjacency_list[s]:
            if (s < t):
                for v in data[s].intersection(data[t]):
                    nn_triangles.append((v, s, t))
                    in_triangles.add(v)
                    in_triangles.add(s)
                    in_triangles.add(t)
                data[t].add(s)
    #set(range(90)).difference(in_triangles)
    return nn_triangles

def confirm_triangles(nn_triangles, list_of_SEACell_assignments_per_cell, count_threshold):
    """
    Check if a nn_triangle is sufficiently represented by the list_of_SEACell_assignments_per_cell. If not, remove them..
    Checking means looking if there are more than count_threshold cells that have non-trivial SEACell assignments to all three SEACells in each nn_triangle.
    """
    counts = np.array([0] * len(nn_triangles))
    confirmed_triangles = []
    removed_triangles = []
    for index, triangle in enumerate(nn_triangles):
        for assignments in list_of_SEACell_assignments_per_cell:
            if(len(set(triangle).difference(set(assignments))) == 0):
                counts[index] += 1
        if(counts[index] > count_threshold):
            confirmed_triangles.append(triangle)
        else:
            removed_triangles.append(triangle)
    removed_seacells = set(range(90)).difference(set(np.array(confirmed_triangles, dtype=object).flatten()))
    return confirmed_triangles, removed_triangles, removed_seacells, counts

def triangles(ad, A):
    """
    Compute the triangles in the graph of SEACells.
    """

    from classes import Data

    # Take the anndata for SEACells
    SEACell_ad = summarise_by_SEACell(ad, A, summarize_layer='raw')
    number_of_seacells = SEACell_ad.shape[0]

    # Compute the list of SEACell assignments per cell
    list_of_SEACell_assignments_per_cell = seacells_by_weights(A, 0.05)
    #print(list_of_SEACell_assignments_per_cell[:10])

    # Normalize SEACells, log transform and compute highly variable genes
    sc.pp.normalize_per_cell(SEACell_ad)
    sc.pp.log1p(SEACell_ad)
    sc.pp.highly_variable_genes(SEACell_ad, n_top_genes=800)

    # Compute principal components and nearest neighbours algorithm and UMAP
    # Here we use 10 components. This number may also be selected by examining variance explained
    sc.tl.pca(SEACell_ad, n_comps=10, use_highly_variable=True)
    sc.pp.neighbors(SEACell_ad, n_neighbors=5)
    sc.tl.umap(SEACell_ad, n_components=2, min_dist=0.5, spread=1.0)

    # Compute the adjacency matrix
    adjacency_matrix = (SEACell_ad.obsp["connectivities"].toarray() != 0).astype(int)
    #print(adjacency_matrix[0])

    # Compute the adjacency list
    list_order = SEACell_ad.obs_names.astype('U').str.split('-').str[1].astype(np.int)
    adjacency_list = np.array([np.where(adjacency_matrix[i] == 1)[0] for i in range(number_of_seacells)], dtype=object)

    # Make sure the SEACells labels are correct and sorted for the adjacency list
    ordered_list = np.empty(number_of_seacells, dtype=object)
    for i in range(number_of_seacells):
        ordered_list[list_order[i]] = np.sort(list_order[adjacency_list[i]])
    adjacency_list = ordered_list
    #print(adjacency_list[0])

    # Compute the nn_triangles
    nn_triangles = compute_nn_triangles(adjacency_list)
    #print(nn_triangles)

    no_of_nn_triangles = int(np.trace(np.linalg.matrix_power(adjacency_matrix, 3)/6))
    assert no_of_nn_triangles == len(nn_triangles), "There was an error in computing the nearest neighbor triangles."

    # Compute the confirmed_triangles
    confirmed_triangles, removed_triangles, removed_seacells, count_matrix = confirm_triangles(nn_triangles, list_of_SEACell_assignments_per_cell, 0)

    list_of_triangles_for_each_seacell = []
    for i in range(number_of_seacells):
        triangles = [triangle for triangle in confirmed_triangles if triangle[0] == i or triangle[1] == i or triangle[2] == i]
        list_of_triangles_for_each_seacell.append(triangles)

    # Right now both confirmed and removed triangles include just numbers. I did not append "SEACell-" to them.
    # Also, each confirmed_triangle is sorted in ascending order. However, the list of confirmed triangles itself is not strictly sorted.

    data = Data(np.array(confirmed_triangles), np.array(removed_triangles), removed_seacells, adjacency_matrix, SEACell_ad.obs_names, SEACell_ad.obsm['X_umap'], adjacency_list, list_of_triangles_for_each_seacell)

    return data


def adjacency_matrix_to_graph(adj_matrix):
    """
    Convert an adjacency matrix to a networkx graph.
    """
    graph = nx.Graph()
    graph.add_nodes_from(range(adj_matrix.shape[0]))
    for i in range(adj_matrix.shape[0]):
        for j in range(i + 1, adj_matrix.shape[0]):
            if adj_matrix[i, j] == 1:
                graph.add_edge(i, j)
    return graph

def annotate_nodes(labels, graph):
    """
    Annotate nodes in a graph with a list of labels.
    """
    for i, label in enumerate(labels):
        graph.nodes[i]['label'] = label
    return graph

def add_coordinates(coords, graph):
    """
    Add coordinates to a graph.
    """
    for i, coord in enumerate(coords):
        graph.nodes[i]['x'] = coord[0]
        graph.nodes[i]['y'] = coord[1]
    return graph


