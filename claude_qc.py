import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scrublet as scr
from anndata import AnnData
from scipy.stats import pearsonr


def qc(adata: AnnData, section: str,
       genes_low=1.0, genes_high=99.5,
       counts_low=1.0, counts_high=99.5,
       min_genes_floor=200, min_counts_floor=500,
       expected_doublet_rate=0.06) -> AnnData:
    """
    Complete QC pipeline matching the paper's methodology.
    """
    print(f'\n{"=" * 60}')
    print(f'Starting QC for AnnData {section}')
    print(f'Initial cells: {adata.n_obs}')
    print(f'{"=" * 60}\n')

    # Step 1: Remove negative probes
    adata = remove_negative_probes(adata, section)

    # Step 2: Remove outlier cells by genes and transcripts
    adata = remove_outlier_cells_by_genes_and_transcripts(
        adata, section,
        genes_low=genes_low, genes_high=genes_high,
        counts_low=counts_low, counts_high=counts_high,
        min_genes_floor=min_genes_floor,
        min_counts_floor=min_counts_floor
    )

    # Step 3: Run Scrublet to remove doublets
    adata = run_scrublet(adata, section, expected_doublet_rate=expected_doublet_rate)

    print(f'\n{"=" * 60}')
    print(f'Finished QC for AnnData {section}')
    print(f'Final cells: {adata.n_obs}')
    print(f'{"=" * 60}\n')

    return adata


def remove_negative_probes(adata: AnnData, section: str) -> AnnData:
    """Remove negative control probes from the dataset."""
    print(f'Checking for negative probes for AnnData {section}')
    neg_mask = adata.var_names.str.startswith("NegPrb")
    n_neg = int(neg_mask.sum())

    if n_neg:
        print(f'Found negative probes: {adata.var_names[neg_mask].tolist()}')
        adata = adata[:, ~neg_mask].copy()
        print(f'Removed {n_neg} negative probes from AnnData {section}')
    else:
        print(f'No negative probes found in AnnData {section}')

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
    """
    Filter cells based on number of genes and transcript counts.
    """
    print(f'\nStarting outlier removal for AnnData {section}')

    # Calculate metrics
    adata.obs["n_counts"] = np.asarray(adata.X.sum(axis=1)).ravel()
    adata.obs["n_genes"] = np.asarray((adata.X > 0).sum(axis=1)).ravel()

    # Determine thresholds
    genes_min, genes_max = choose_bounds(
        adata.obs["n_genes"], genes_low, genes_high, floor=min_genes_floor
    )
    counts_min, counts_max = choose_bounds(
        adata.obs["n_counts"], counts_low, counts_high, floor=min_counts_floor
    )

    print(f'Thresholds - Genes: [{genes_min:.0f}, {genes_max:.0f}], '
          f'Counts: [{counts_min:.0f}, {counts_max:.0f}]')

    # Create visualization
    fig, ax = plt.subplots(1, 2, figsize=(12, 4))

    ax[0].hist(adata.obs["n_counts"], bins=100, alpha=0.7, edgecolor='black')
    ax[0].axvline(counts_min, ls="--", color='red', linewidth=2, label=f'Min: {counts_min:.0f}')
    ax[0].axvline(counts_max, ls="--", color='red', linewidth=2, label=f'Max: {counts_max:.0f}')
    ax[0].set_xlabel("Total counts per cell")
    ax[0].set_ylabel("Number of cells")
    ax[0].set_title(f"{section} — Transcript counts distribution")
    ax[0].legend()

    ax[1].hist(adata.obs["n_genes"], bins=100, alpha=0.7, color="orange", edgecolor='black')
    ax[1].axvline(genes_min, ls="--", color='red', linewidth=2, label=f'Min: {genes_min:.0f}')
    ax[1].axvline(genes_max, ls="--", color='red', linewidth=2, label=f'Max: {genes_max:.0f}')
    ax[1].set_xlabel("Detected genes per cell")
    ax[1].set_ylabel("Number of cells")
    ax[1].set_title(f"{section} — Gene detection distribution")
    ax[1].legend()

    plt.tight_layout()
    plt.show()

    # Apply filters
    filter_condition = (
            (adata.obs["n_genes"] >= genes_min) & (adata.obs["n_genes"] <= genes_max) &
            (adata.obs["n_counts"] >= counts_min) & (adata.obs["n_counts"] <= counts_max)
    )

    n_filtered = (~filter_condition).sum()
    adata_filtered = adata[filter_condition].copy()

    print(f"Removed {n_filtered} outlier cells ({n_filtered / adata.n_obs:.2%})")
    print(f"Retained {adata_filtered.n_obs} / {adata.n_obs} cells")

    return adata_filtered


def choose_bounds(x, low=1.0, high=99.5, floor=None, ceil=None):
    """Calculate bounds based on percentiles with optional floor/ceiling."""
    low_val, high_val = np.percentile(x, [low, high])
    if floor is not None:
        low_val = max(low_val, floor)
    if ceil is not None:
        high_val = min(high_val, ceil)
    return float(low_val), float(high_val)


def run_scrublet(adata: AnnData, section: str,
                 expected_doublet_rate=0.06,
                 manual_threshold=None) -> AnnData:
    """
    Run Scrublet to identify and remove doublets.
    """
    print(f'\nRunning Scrublet for AnnData {section}')
    print(f'Expected doublet rate: {expected_doublet_rate:.1%}')

    counts_matrix = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X

    scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=expected_doublet_rate)
    doublet_scores, predicted_doublets = scrub.scrub_doublets()

    adata.obs["doublet_score"] = doublet_scores
    adata.obs["predicted_doublet"] = predicted_doublets

    # Visualization
    fig, ax = plt.subplots(1, 1, figsize=(8, 4))
    ax.hist(doublet_scores, bins=100, color="steelblue", alpha=0.7, edgecolor='black')

    if manual_threshold is not None:
        ax.axvline(manual_threshold, ls="--", color='red', linewidth=2,
                   label=f'Manual threshold: {manual_threshold:.3f}')
        adata.obs["predicted_doublet"] = adata.obs["doublet_score"] > manual_threshold
    else:
        threshold = scrub.threshold_
        ax.axvline(threshold, ls="--", color='red', linewidth=2,
                   label=f'Auto threshold: {threshold:.3f}')

    ax.set_xlabel("Doublet score")
    ax.set_ylabel("Number of cells")
    ax.set_title(f"{section} — Scrublet doublet score distribution")
    ax.legend()
    plt.tight_layout()
    plt.show()

    # Remove doublets
    n_before = adata.n_obs
    adata_clean = adata[~adata.obs["predicted_doublet"]].copy()
    n_after = adata_clean.n_obs
    n_removed = n_before - n_after

    print(f"Removed {n_removed} doublets ({n_removed / n_before:.2%})")

    return adata_clean


def find_correlated_plaque(adata: AnnData, section: str,
                           plaque_col='Mean.BetaAmyloid',
                           plaque_threshold=0.0,
                           return_correlations=True):
    """
    Correlate gene expression with plaque overlap.
    Should be run AFTER adding negative probes back for analysis.

    Args:
        adata: AnnData object (with negative probes included)
        section: Section name for display
        plaque_col: Column name for plaque measurements
        plaque_threshold: Minimum plaque overlap to include cells
        return_correlations: Whether to return correlation dataframe
    """
    print(f'\n{"=" * 60}')
    print(f'Analyzing gene-plaque correlations for {section}')
    print(f'{"=" * 60}\n')

    if plaque_col not in adata.obs.columns:
        print(f"WARNING: Column '{plaque_col}' not found in adata.obs")
        print(f"Available columns: {adata.obs.columns.tolist()}")
        return None

    # Filter for cells overlapping with plaques
    plaque_mask = adata.obs[plaque_col] > plaque_threshold
    n_plaque_cells = plaque_mask.sum()

    print(f"Cells overlapping with plaques (>{plaque_threshold}): {n_plaque_cells} / {adata.n_obs}")

    if n_plaque_cells == 0:
        print("No cells overlap with plaques. Cannot compute correlations.")
        return None

    adata_plaque = adata[plaque_mask].copy()

    # Calculate correlations
    correlations = []
    for gene in adata_plaque.var_names:
        gene_counts = adata_plaque[:, gene].X.toarray().ravel()
        plaque_values = adata_plaque.obs[plaque_col].values

        if np.std(gene_counts) > 0 and np.std(plaque_values) > 0:
            r, p_val = pearsonr(gene_counts, plaque_values)
            correlations.append({
                'gene': gene,
                'R': r,
                'R2': r ** 2,
                'p_value': p_val
            })
        else:
            correlations.append({
                'gene': gene,
                'R': 0,
                'R2': 0,
                'p_value': 1.0
            })

    correlation_df = pd.DataFrame(correlations)
    correlation_df = correlation_df.sort_values("R2", ascending=False)

    # Display top correlations
    print(f"\nTop 20 genes/probes correlated with plaque overlap:")
    print(correlation_df.head(20).to_string(index=False))

    # Identify negative probes in top correlations
    neg_probes = correlation_df[correlation_df['gene'].str.startswith('NegPrb')]
    if len(neg_probes) > 0:
        print(f"\nNegative probes found:")
        print(neg_probes.to_string(index=False))

        # Identify genes to potentially remove (higher R2 than top negative probes)
        if len(neg_probes) >= 2:
            top_neg_r2 = neg_probes.iloc[1]['R2']  # Second highest
            genes_to_remove = correlation_df[
                (correlation_df['R2'] > top_neg_r2) &
                (~correlation_df['gene'].str.startswith('NegPrb'))
                ]
            print(f"\nGenes with R² > {top_neg_r2:.3f} (top 2 neg probes threshold): {len(genes_to_remove)}")
            print(genes_to_remove[['gene', 'R2']].to_string(index=False))

    # Visualization
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))

    # Separate negative probes and genes
    is_neg = correlation_df['gene'].str.startswith('NegPrb')

    ax.scatter(range(len(correlation_df[~is_neg])),
               correlation_df[~is_neg]['R2'],
               alpha=0.5, s=20, label='Genes', color='steelblue')

    neg_indices = correlation_df[is_neg].index
    ax.scatter([correlation_df.index.get_loc(i) for i in neg_indices],
               correlation_df.loc[neg_indices, 'R2'],
               alpha=0.8, s=50, label='Negative probes', color='red', marker='X')

    ax.set_xlabel("Genes (ranked by R²)")
    ax.set_ylabel("R² (correlation with plaque overlap)")
    ax.set_title(f"{section} — Gene expression correlation with plaque overlap")
    ax.legend()
    plt.tight_layout()
    plt.show()

    if return_correlations:
        return correlation_df

    return None


def remove_genes_by_plaque_correlation(adata: AnnData,
                                       correlation_df: pd.DataFrame,
                                       r2_threshold: float) -> AnnData:
    """
    Remove genes that have higher plaque correlation than threshold.

    Args:
        adata: AnnData object
        correlation_df: DataFrame from find_correlated_plaque
        r2_threshold: R² threshold (e.g., from top 2 negative probes)
    """
    genes_to_remove = correlation_df[
        (correlation_df['R2'] > r2_threshold) &
        (~correlation_df['gene'].str.startswith('NegPrb'))
        ]['gene'].tolist()

    print(f"\nRemoving {len(genes_to_remove)} genes with R² > {r2_threshold:.3f}")

    genes_to_keep = [g for g in adata.var_names if g not in genes_to_remove]
    adata_filtered = adata[:, genes_to_keep].copy()

    print(f"Genes remaining: {adata_filtered.n_vars} / {adata.n_vars}")

    return adata_filtered