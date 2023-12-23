import numpy as np
import scanpy as sc
from tqdm import tqdm
import networkx as nx
from math import sqrt
from numpy import random
from numpy import linalg
from seacell_computer import sparsify_assignments
from seacell_computer import summarise_by_SEACell

def triangle_info(ad, A, seacell_1 = None, seacell_2 = None, seacell_3 = None):
    """
    Given three SEACells, return the adjacency matrix, labels, and coordinates, and strength list for the triangle that they form with the cells inside.
    :param ad: (AnnData) anndata object containing single cell data
    :param A: (csr_matrix) of shape n_cells x n_SEACells containing assignment weights
    :param seacell_1: (int) index of first SEACell
    :param seacell_2: (int) index of second SEACell
    :param seacell_3: (int) index of third SEACell
    :return: (np.array) of the adjacency matrix of the triangle network, (list) of labels of the triangle, (list) of coordinates of the triangle, (list) of strength of the triangle
    """
    v1 = [0,sqrt(3)]
    v2 = [1,0]
    v3 = [-1,0]
    tri_coords = [v1,v2,v3]
    strength = [1, 1, 1]

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
        
        #if(linalg.norm([coordinate[0] - v1[0],coordinate[1] - v1[1]]) < 0.0001 or linalg.norm([coordinate[0] - v2[0],coordinate[1] - v2[1]]) < 0.0001 or linalg.norm([coordinate[0] - v3[0],coordinate[1] - v3[1]]) < 0.0001):
        if(coordinate == v1 or coordinate == v2 or coordinate == v3):
            coordinate = [coordinate[0] + 0.005 * random.rand() - 0.0025, coordinate[1] +  0.005 * random.rand() - 0.0025]
        
        tri_coords.append(coordinate)
        strength.append(int(sum/0.05)/20)
        tri_labels.append(ad.obs.index[i])

    tri_adj_matrix = np.zeros(shape=(len(tri_coords), len(tri_coords)), dtype=np.int8)
    tri_adj_matrix[0, 1:3], tri_adj_matrix[1:3, 0], tri_adj_matrix[2, 1], tri_adj_matrix[1, 2] = 1, 1, 1, 1
    return tri_adj_matrix, tri_labels, tri_coords, strength

def seacells_by_weights(A, threshold):
    """
    Return a list of proper SEACells assignments for each cell by thresholding the weights.
    :param A: (csr_matrix) of shape n_cells x n_SEACells containing assignment weights
    :param threshold: (float) threshold above which the assignment weight is considered non-trivial
    :return: (list) of SEACell assignments for each cell
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
    Compute the triangles in the SEACElls graph using the adjacency list that the nearest neighbor method provides.
    :param adjacency_list: (list) of lists of neighbors for each SEACell
    :return: (list) of triangles in the SEACell graph
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
    Returns the list of triangles that are sufficiently represented by the SEACell assignments.
    :param nn_triangles: (list) of triangles in the SEACell graph computed by the nearest neighbor method
    :param list_of_SEACell_assignments_per_cell: (list) of SEACell assignments for each cell
    :param count_threshold: (int) threshold for the number of cells assigned to a SEACell triangle, above which that triangle is considered sufficiently represented
    :return: (list) of confirmed triangles in the SEACell graph
    """
    counts = np.array([0] * len(nn_triangles))
    confirmed_triangles = []
    removed_triangles = []
    for index, triangle in enumerate(nn_triangles):
        for assignments in list_of_SEACell_assignments_per_cell:
            if(len(set(triangle).difference(set(assignments))) == 0):
                counts[index] += 1
        if(counts[index] >= 2):
            confirmed_triangles.append(triangle)
        else:
            removed_triangles.append(triangle)
    removed_seacells = set(range(90)).difference(set(np.array(confirmed_triangles, dtype=object).flatten()))
    return confirmed_triangles, removed_triangles, removed_seacells, counts

def triangles(ad, A, SEACell_ad):
    """
    Compute the triangles in the graph of SEACells.
    :param ad: (AnnData) anndata object containing single cell data
    :param A: (csr_matrix) of shape n_cells x n_SEACells containing assignment weights
    :return: (Data) object containing the adjacency matrix, labels, coordinates of the graph of SEACells, the list of confirmed triangles of SEACells, and the other information in Data object
    """
    
    from classes import Data


    ''' # Take the anndata for SEACells
    #SEACell_ad = summarise_by_SEACell(ad, A, summarize_layer='raw')
    if (seacell_ad is None):
        SEACell_ad = sc.read("data/seacell_anndata.h5ad")
    else:
        SEACell_ad = seacell_ad'''
        
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
    list_order = SEACell_ad.obs_names.astype('U').str.split('-').str[1].astype(int)
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
    print('HEre are the triangle numbers:')
    print(no_of_nn_triangles)
    print(len(nn_triangles))

    # assert no_of_nn_triangles == len(nn_triangles), "There was an error in computing the nearest neighbor triangles."

    # Compute the confirmed_triangles
    confirmed_triangles, removed_triangles, removed_seacells, count_matrix = confirm_triangles(nn_triangles, list_of_SEACell_assignments_per_cell, 3)

    print('Here are the counts:')
    print(count_matrix)
    print('confirmed triangles, removed triangles, removed seacells')
    print(len(confirmed_triangles))
    print(len(removed_triangles))
    print(len(removed_seacells))

    list_of_triangles_for_each_seacell = []
    for i in range(number_of_seacells):
        triangles = [triangle for triangle in confirmed_triangles if triangle[0] == i or triangle[1] == i or triangle[2] == i]
        list_of_triangles_for_each_seacell.append(triangles)

    print('Here are the triangles for each seacell:')
    print(list_of_triangles_for_each_seacell)

    # Right now both confirmed and removed triangles include just numbers. I did not append "SEACell-" to them.
    # Also, each confirmed_triangle is sorted in ascending order. However, the list of confirmed triangles itself is not strictly sorted.

    data = Data(np.array(confirmed_triangles), np.array(removed_triangles), removed_seacells, adjacency_matrix, SEACell_ad.obs_names, SEACell_ad.obsm['X_umap'], adjacency_list, list_of_triangles_for_each_seacell)

    return data




def adjacency_matrix_to_graph(adj_matrix):
    """
    Convert an adjacency matrix to a networkx graph.
    :param adj_matrix: (np.array) of shape n_nodes x n_nodes containing the adjacency matrix
    :return: (networkx.Graph) graph corresponding to the adjacency matrix
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
    :param labels: (list) of labels for each node
    :param graph: (networkx.Graph) graph to annotate
    :return: (networkx.Graph) annotated graph
    """
    for i, label in enumerate(labels):
        graph.nodes[i]['label'] = label
    return graph

def add_coordinates(coords, graph):
    """
    Add coordinates to a graph.
    :param coords: (list) of coordinates for each node
    :param graph: (networkx.Graph) graph to annotate
    :return: (networkx.Graph) graph with coordinates
    """
    for i, coord in enumerate(coords):
        graph.nodes[i]['x'] = coord[0]
        graph.nodes[i]['y'] = coord[1]
    return graph


