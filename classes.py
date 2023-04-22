class Data:
    def __init__(self, confirmed_triangles, removed_triangles, removed_seacells, adjacency_matrix, seacell_labels, UMAP_coords, adjacency_list, list_of_triangles_for_each_seacell):
        self.confirmed_triangles = confirmed_triangles
        self.removed_triangles = removed_triangles
        self.removed_seacells = removed_seacells
        self.adjacency_matrix = adjacency_matrix
        self.seacell_labels = seacell_labels
        self.UMAP_coords = UMAP_coords
        self.adjacency_list = adjacency_list
        self.list_of_triangles_for_each_seacell = list_of_triangles_for_each_seacell


class Input:
    def __init__(self, ad, model):
        self.ad = ad
        self.model = model
