import anndata as ad
import squidpy as sq
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from qc import qc, find_correlated_plaque

DATA_PATH = "D:/thesis_research/resources"


def load_adata(section_id: str) -> ad.AnnData:
    return sq.read.nanostring(
    path=DATA_PATH,
    counts_file=f'{section_id}_exprMat_file.csv',
    meta_file=f'{section_id}_metadata_file.csv',
    fov_file=f'{section_id}_fov_positions_file.csv')


def run():
    adata_s10 = load_adata(section_id='GSM8199188_ID61-ID62_S10')
    adata_s18 = load_adata(section_id='GSM8199189_ID67-ID68_S18')
    find_correlated_plaque(adata_s10)
    find_correlated_plaque(adata_s18)

    qc(adata = adata_s10, section='s10')
    qc(adata = adata_s18, section='s18')


run()















#

#adata_s10.obs['section'] = 'S10'
#adata_s18.obs['section'] = 'S18'
# adata_all = adata_s10.concatenate(adata_s18, batch_key="section", index_unique=None)

# print(adata_all)
# print(adata_all.obs.head())
# print(adata_all.var.head())
#
# print(f"Cells: {adata_all.n_obs:,}")
# print(f"Genes: {adata_all.n_vars:,}")
# print("Negative probes present:",
#       [g for g in adata_all.var_names if 'NegPrb' in g][:10])

# Calculate QC metrics



# fig, axs = plt.subplots(1, 2, figsize=(10,4))
# axs[0].hist(adata_all.obs['n_counts'], bins=100, color='steelblue')
# axs[0].set_xlabel('Total transcripts per cell')
# axs[1].hist(adata_all.obs['n_genes'], bins=100, color='coral')
# axs[1].set_xlabel('Genes detected per cell')
# plt.show()