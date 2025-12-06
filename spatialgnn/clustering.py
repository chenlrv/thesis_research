from anndata import AnnData


def cluster_from_gnn(
        adata: AnnData,
        gnn_key: str = "X_gnn_dgi",
        n_neighbors: int = 15,
        resolution: float = 0.6,
        plot: bool = True,
):
    """
    Cluster ALL cells using GNN embeddings.
    """
    import scanpy as sc

    if gnn_key not in adata.obsm:
        raise KeyError(f"{gnn_key} not found in adata.obsm.")

    adata.obsm["X_gnn_used"] = adata.obsm[gnn_key]

    sc.pp.neighbors(adata, use_rep="X_gnn_used", n_neighbors=n_neighbors)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=resolution, key_added="gnn_cluster")

    if plot:
        sc.pl.umap(adata, color="gnn_cluster")

    return adata
