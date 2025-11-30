import scanpy as sc
from anndata import AnnData
import numpy as np
import scipy.sparse as sp

def clustering(adata: AnnData) -> None:
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    sc.pp.scale(adata, max_value=10)

    sc.tl.pca(adata, n_comps=30)

    sc.pp.neighbors(adata, n_pcs=20)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=0.45, key_added="cluster")

    sc.pl.umap(adata, color=["cluster"], legend_loc="on data")

    _label_clusters(adata)

    _filter_microglia(adata)


def _label_clusters(adata: AnnData):
    markers = _get_marker_genes()

    for cell_type, genes in markers.items():
        sc.tl.score_genes(adata, gene_list=genes, score_name=f"{cell_type}_score")

    cluster_labels = {}
    for cluster in sorted(adata.obs["cluster"].unique(), key=int):
        sub = adata[adata.obs["cluster"] == cluster]
        mean_scores = sub.obs[[f"{cell_type}_score" for cell_type in markers]].mean()
        cluster_labels[cluster] = mean_scores.idxmax().replace("_score", "")

    adata.obs["celltype"] = adata.obs["cluster"].map(cluster_labels)

    sc.pl.umap(adata, color="celltype", legend_loc="on data")


def _filter_microglia(adata: AnnData):
    micro = adata[adata.obs["celltype"] == "Microglia"].copy()
    print("Microglia shape before QC:", micro.shape)

    print("micro shape before PCA:", micro.shape)

    neg_probes = [g for g in adata.var_names
                  if "neg" in g.lower()]

    sc.tl.pca(micro, n_comps=30)
    sc.pp.neighbors(micro, n_pcs=15)
    sc.tl.umap(micro)
    sc.tl.leiden(micro, resolution=0.4, key_added="micro_cluster")

    sc.pl.umap(micro, color="micro_cluster", legend_loc="on data")

    df = micro.to_df()

    cluster_means = (
        df[neg_probes]
        .groupby(micro.obs["micro_cluster"])
        .mean()
    )

    bad_cluster = cluster_means.mean(axis=1).idxmax()

    bad_cells = micro.obs[micro.obs["micro_cluster"] == bad_cluster].index
    len(bad_cells)

    adata_clean = adata[~adata.obs.index.isin(bad_cells)].copy()
    adata_clean


def _get_marker_genes() -> dict[str, list[str]]:
    return {
        "ChoroidPlexus": ["Ttr"],
        "Ependymal": ["Foxj1"],
        "Oligodendrocytes": ["Mbp", "Mog", "Mag", "Plp1", "Bcas1", "Ugt8a", "Myrf", "Olig1", "Olig2"],
        "Astrocytes": ["Gfap", "Aqp4", "Aldh1l1", "Slc1a3", "S100b", "Igfbp7"],
        "Microglia": ["Cx3cr1", "P2ry12", "Trem2", "C1qa", "C1qb", "C1qc", "Csf1r", "Tmem119", "Sall1"],
        "OPC": ["Pdgfra", "Olig2"],
        "Endothelia": ["Cldn5", "Pecam1", "Flt1", "Kdr", "Esam", "Emcn", "Vegfa"],

        "Pericytes": ["Pdgfrb", "Rgs5", "Acta2", "Tagln", "Col1a1", "Col1a2", "Des"],
        "GABAergic": ["Gad1", "Gad2", "Slc32a1", "Pvalb", "Sst", "Vip", "Crh"],
        "Glutamatergic": ["Slc17a7", "Rbfox3", "Camk2a", "Thy1"],
    }
