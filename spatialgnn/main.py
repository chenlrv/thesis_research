# core
import numpy as np
import scanpy as sc
from anndata import AnnData

# graph & DL
import torch
from torch import nn
import torch.nn.functional as F
from torch_geometric.data import Data
from torch_geometric.nn import GCNConv
from torch_geometric.nn.models import DGI

# kNN graph
from sklearn.neighbors import NearestNeighbors

def prepare_node_features(adata: AnnData, n_hvg: int = 2000, n_pcs: int = 50) -> np.ndarray:
    """
    Take AnnData, select HVGs, run PCA, return node features (cells × n_pcs).
    Stores results in adata so we can reuse.
    """
    # Work on a copy to avoid messing with your original unintentionally
    adata_tmp = adata.copy()

    # Highly variable genes
    sc.pp.highly_variable_genes(adata_tmp, n_top_genes=n_hvg, subset=True)

    # Normalize & log
    sc.pp.normalize_total(adata_tmp, target_sum=1e4)
    sc.pp.log1p(adata_tmp)

    # PCA
    sc.pp.scale(adata_tmp, max_value=10)
    sc.tl.pca(adata_tmp, n_comps=n_pcs)

    # Store in original adata
    adata.obsm["X_pca_gnn"] = adata_tmp.obsm["X_pca"]

    return adata.obsm["X_pca_gnn"]


def build_spatial_knn_graph(adata: AnnData, k: int = 8,
                            x_key: str = "CenterX_global_px",
                            y_key: str = "CenterY_global_px"):
    """
    Build a kNN graph based on spatial coordinates.
    Returns edge_index (2 × num_edges) as torch.LongTensor.
    """
    coords = np.vstack([adata.obs[x_key].values, adata.obs[y_key].values]).T  # (n_cells, 2)
    n_cells = coords.shape[0]

    nn = NearestNeighbors(n_neighbors=k + 1, metric="euclidean")  # +1 for self
    nn.fit(coords)
    dists, indices = nn.kneighbors(coords)  # indices: (n_cells, k+1)

    # Build edge list (i -> j)
    # skip index 0 because it's the cell itself
    edge_src = []
    edge_dst = []
    for i in range(n_cells):
        neighbors = indices[i, 1:]  # ignore self
        for j in neighbors:
            edge_src.append(i)
            edge_dst.append(j)

    edge_index = np.vstack([edge_src, edge_dst])  # shape (2, num_edges)

    # Make it undirected by adding reverse edges
    edge_index_rev = edge_index[::-1, :]
    edge_index = np.hstack([edge_index, edge_index_rev])

    # Optionally, remove duplicates
    edge_index = torch.as_tensor(edge_index, dtype=torch.long)
    edge_index = torch.unique(edge_index, dim=1)  # unique columns

    return edge_index


def adata_to_pyg_data(adata: AnnData,
                      feature_key: str = "X_pca_gnn",
                      x_key: str = "CenterX_global_px",
                      y_key: str = "CenterY_global_px",
                      k: int = 8) -> Data:
    """
    Convert AnnData to PyG Data object:
    - x = node features (PCA)
    - edge_index = spatial kNN graph
    - pos = spatial coordinates (for potential spatial encodings)
    """
    if feature_key not in adata.obsm:
        raise ValueError(f"{feature_key} not found in adata.obsm. Run prepare_node_features first.")

    x = torch.as_tensor(adata.obsm[feature_key], dtype=torch.float32)
    edge_index = build_spatial_knn_graph(adata, k=k, x_key=x_key, y_key=y_key)

    pos = torch.as_tensor(
        np.vstack([adata.obs[x_key].values,
                   adata.obs[y_key].values]).T,
        dtype=torch.float32
    )

    data = Data(x=x, edge_index=edge_index, pos=pos)
    return data


X_pca = prepare_node_features(adata_all, n_hvg=2000, n_pcs=50)
print(X_pca.shape)  # (n_cells, 50)
edge_index = build_spatial_knn_graph(adata_all, k=8,
                                     x_key="CenterX_global_px",
                                     y_key="CenterY_global_px")
print(edge_index.shape)  # (2, num_edges)