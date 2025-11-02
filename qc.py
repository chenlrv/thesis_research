import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scrublet as scr
from anndata import AnnData
from scipy.stats import pearsonr

def qc(adata: AnnData, section: str) -> AnnData:
    print(f'Starting QC for AnnData {section}')
    adata = remove_negative_probes(adata, section)
    adata = remove_outlier_cells_by_genes_and_transcripts(adata, section)
    adata = run_scrublet(adata, section)
    print(f'Finished QC for AnnData {section}')
    return adata


def remove_negative_probes(adata: AnnData, section: str) -> AnnData:
    print(f'Checking for negative probes for Anndata {section}')
    neg_mask = adata.var_names.str.startswith("NegPrb")
    n_neg = int(neg_mask.sum())
    if n_neg:
        adata = adata[:, ~neg_mask].copy()
    print(f'Removed {n_neg} negative probes from Anndata {section}')
    return adata


def remove_outlier_cells_by_genes_and_transcripts(
        adata: AnnData,
        section: str,
        genes_low=1.0,
        genes_high=99.5,
        counts_low=1.0,
        counts_high=99.5,
        min_genes_floor=200,
        min_counts_floor=500) -> AnnData:
    print(f'Starting outlier removal for Anndata {section}')

    adata.obs["n_counts"] = np.asarray(adata.X.sum(axis=1)).ravel()
    adata.obs["n_genes"] = np.asarray((adata.X > 0).sum(axis=1)).ravel()

    genes_min, genes_max = choose_bounds(adata.obs["n_genes"], genes_low, genes_high, floor=min_genes_floor)
    counts_min, counts_max = choose_bounds(adata.obs["n_counts"], counts_low, counts_high, floor=min_counts_floor)

    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    ax[0].hist(adata.obs["n_counts"], bins=100)
    ax[0].axvline(counts_min, ls="--")
    ax[0].axvline(counts_max, ls="--")
    ax[0].set_title(f"{section} — n_counts")
    ax[1].hist(adata.obs["n_genes"], bins=100, color="orange")
    ax[1].axvline(genes_min, ls="--")
    ax[1].axvline(genes_max, ls="--")
    ax[1].set_title(f"{section} — n_genes")
    plt.show()

    filter_condition = (
            (adata.obs["n_genes"] >= genes_min) & (adata.obs["n_genes"] <= genes_max) &
            (adata.obs["n_counts"] >= counts_min) & (adata.obs["n_counts"] <= counts_max)
    )
    adata_filtered = adata[filter_condition].copy()
    print(f"AnnData {section}: kept {adata_filtered.n_obs} / {adata.n_obs} "
          f"(genes {genes_min:.0f}-{genes_max:.0f}, counts {counts_min:.0f}-{counts_max:.0f})")
    thresholds = dict(gmin=genes_min, gmax=genes_max, cmin=counts_min, cmax=counts_max)
    return adata_filtered


def choose_bounds(x, low=1.0, high=99.5, floor=None, ceil=None):
    low, high = np.percentile(x, [low, high])
    if floor is not None: low = max(low, floor)
    if ceil is not None: high = min(high, ceil)
    return float(low), float(high)


def run_scrublet(adata: AnnData, section: str, expected_doublet_rate=0.06, manual_threshold=None) -> AnnData:
    """ Run Scrublet to identify doublets in the AnnData object."""
    print(f'Running Scrublet for Anndata {section}')

    counts_matrix = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X

    scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=expected_doublet_rate)
    doublet_scores, predicted_doublets = scrub.scrub_doublets()

    adata.obs["doublet_score"] = doublet_scores
    adata.obs["predicted_doublet"] = predicted_doublets

    plt.figure(figsize=(5, 3))
    plt.hist(doublet_scores, bins=100, color="steelblue")
    plt.title(f"{section} Scrublet doublet scores")
    plt.xlabel("Doublet score")
    plt.ylabel("Cell count")
    plt.show()

    if manual_threshold is not None:
        adata.obs["predicted_doublet"] = adata.obs["doublet_score"] > manual_threshold
        print(f"Manual threshold set at {manual_threshold:.3f}")

    n_before = adata.n_obs
    adata = adata[~adata.obs["predicted_doublet"]].copy()
    n_after = adata.n_obs

    print(f"Removed from AnnData {section}: {n_before - n_after} doublets ({(n_before - n_after) / n_before:.2%})")
    return adata


def find_correlated_plaque(adata: AnnData):
    plaque_col = 'Mean.BetaAmyloid'

    correlations = []
    for gene in adata.var_names:
        gene_counts = adata[:, gene].X.toarray().ravel()
        mean_plaque_values = adata.obs[plaque_col].values
        if np.std(gene_counts) > 0 and np.std(mean_plaque_values) > 0:
            r, _ = pearsonr(gene_counts, mean_plaque_values)
            correlations.append(r ** 2)
        else:
            correlations.append(0)

    correlation = pd.DataFrame({"gene": adata.var_names, "R2": correlations})
    correlation.sort_values("R2", ascending=False).head(10)



# --- run per section (after dropping cell_ID==0 earlier) ---
# s10_filt, s10_bounds = qc_per_section(adata_s10, "S10")
# s18_filt, s18_bounds = qc_per_section(adata_s18, "S18")
#
# # re-combine
# adata_all_qc1 = s10_filt.concatenate(s18_filt, batch_key="section", batch_categories=["S10","S18"], index_unique=None)
# print("After basic QC:", adata_all_qc1.n_obs)
