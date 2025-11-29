def cluster_macrophages_from_gnn(adata: AnnData,
                                 gnn_key: str = "X_gnn_dgi",
                                 celltype_col: str = "cell_type",
                                 macrophage_label: str = "Macrophage",
                                 resolution: float = 0.6):
    """
    Subset macrophages, cluster them in GNN embedding space.
    """
    # subset
    mg = adata[adata.obs[celltype_col] == macrophage_label].copy()

    if gnn_key not in mg.obsm:
        raise ValueError(f"{gnn_key} not found in obsm. Did you run get_gnn_embeddings?")

    mg.obsm["X_gnn_used"] = mg.obsm[gnn_key]

    # neighbors & clustering
    sc.pp.neighbors(mg, use_rep="X_gnn_used", n_neighbors=15)
    sc.tl.umap(mg)
    sc.tl.leiden(mg, resolution=resolution, key_added="mg_gnn_cluster")

    return mg