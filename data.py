class Data:
    def __init__(self, validated_triangles, removed_triangles, removed_seacells, adjacency_matrix, celltypes, UMAP_coords):
        self.validated_triangles = validated_triangles
        self.removed_triangles = removed_triangles
        self.removed_seacells = removed_seacells
        self.adjacency_matrix = adjacency_matrix
        self.celltypes = celltypes
        self.UMAP_coords = UMAP_coords

class Input:
    def __init__(self, ad, model):
        self.ad = ad
        self.model = model