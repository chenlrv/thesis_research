import anndata as ad
import squidpy as sq

from qc import qc, find_correlated_plaque


def load_adata(section_id: str, data_path: str) -> ad.AnnData:
    """Load CosMx data for a single section."""
    return sq.read.nanostring(
        path=data_path,
        counts_file=f'{section_id}_exprMat_file.csv',
        meta_file=f'{section_id}_metadata_file.csv',
        fov_file=f'{section_id}_fov_positions_file.csv'
    )


def complete_qc_workflow(data_path: str):
    """
    Complete QC workflow matching the paper's methodology.

    Steps:
    1. Load data for both sections
    2. Remove cells with cell_ID == 0 (if present)
    3. Basic QC (remove negative probes, outliers, doublets) - per section
    4. Analyze plaque correlation (WITH negative probes added back)
    5. Remove genes with high plaque correlation
    6. Additional microglial QC (if needed)
    """

    print("=" * 80)
    print("STEP 1: Loading data")
    print("=" * 80)

    # Load both sections
    adata_s10 = load_adata('GSM8199188_ID61-ID62_S10', data_path)
    adata_s18 = load_adata('GSM8199189_ID67-ID68_S18', data_path)

    print(f"\nS10 initial: {adata_s10.n_obs} cells, {adata_s10.n_vars} genes")
    print(f"S18 initial: {adata_s18.n_obs} cells, {adata_s18.n_vars} genes")
    print(f"Total initial: {adata_s10.n_obs + adata_s18.n_obs} cells")

    # Check for cell_ID == 0
    for name, adata in [('S10', adata_s10), ('S18', adata_s18)]:
        if 'cell_ID' in adata.obs.columns:
            n_zero = (adata.obs['cell_ID'] == 0).sum()
            if n_zero > 0:
                print(f"\nWARNING: {name} has {n_zero} cells with cell_ID == 0")
                print("Consider removing: adata = adata[adata.obs['cell_ID'] != 0].copy()")

    print("\n" + "=" * 80)
    print("STEP 2: Basic QC (per section)")
    print("=" * 80)

    # Save negative probes for later correlation analysis
    neg_mask_s10 = adata_s10.var_names.str.startswith("NegPrb")
    neg_mask_s18 = adata_s18.var_names.str.startswith("NegPrb")

    neg_probes_s10 = adata_s10[:, neg_mask_s10].copy()
    neg_probes_s18 = adata_s18[:, neg_mask_s18].copy()

    # Run basic QC (removes negative probes, outliers, doublets)
    # Adjust these parameters to match your data
    adata_s10_qc = qc(
        adata_s10,
        section='S10',
        genes_low=2.0,  # Try different percentiles
        genes_high=98.0,
        counts_low=2.0,
        counts_high=98.0,
        min_genes_floor=200,
        min_counts_floor=500,
        expected_doublet_rate=0.06
    )

    adata_s18_qc = qc(
        adata_s18,
        section='S18',
        genes_low=2.0,
        genes_high=98.0,
        counts_low=2.0,
        counts_high=98.0,
        min_genes_floor=200,
        min_counts_floor=500,
        expected_doublet_rate=0.06
    )

    print(f"\nAfter basic QC:")
    print(f"S10: {adata_s10_qc.n_obs} cells")
    print(f"S18: {adata_s18_qc.n_obs} cells")
    print(f"Total: {adata_s10_qc.n_obs + adata_s18_qc.n_obs} cells")

    print("\n" + "=" * 80)
    print("STEP 3: Plaque correlation analysis")
    print("=" * 80)

    # Add negative probes back for correlation analysis
    # Match cells that passed QC
    neg_probes_s10_qc = neg_probes_s10[adata_s10_qc.obs_names].copy()
    neg_probes_s18_qc = neg_probes_s18[adata_s18_qc.obs_names].copy()

    # Concatenate with QC'd data
    adata_s10_with_neg = ad.concat([adata_s10_qc, neg_probes_s10_qc], axis=1)
    adata_s18_with_neg = ad.concat([adata_s18_qc, neg_probes_s18_qc], axis=1)

    # Run correlation analysis
    # Note: Check what the plaque column is actually called in your metadata
    corr_s10 = find_correlated_plaque(
        adata_s10_with_neg,
        section='S10',
        plaque_col='Mean.BetaAmyloid',  # Adjust this name if needed
        plaque_threshold=0.0
    )

    corr_s18 = find_correlated_plaque(
        adata_s18_with_neg,
        section='S18',
        plaque_col='Mean.BetaAmyloid',
        plaque_threshold=0.0
    )

    print("\n" + "=" * 80)
    print("STEP 4: Remove genes with high plaque correlation")
    print("=" * 80)

    # Get R² threshold from top 2 negative probes
    # You'll need to manually check the correlation results and set this
    # Based on paper: NegPrb1 R²=0.65, NegPrb7 R²=0.47
    # They used the second highest (0.47) as threshold

    if corr_s10 is not None:
        neg_probes_df = corr_s10[corr_s10['gene'].str.startswith('NegPrb')]
        if len(neg_probes_df) >= 2:
            r2_threshold = neg_probes_df.iloc[1]['R2']
            print(f"\nS10 R² threshold (from 2nd highest neg probe): {r2_threshold:.3f}")
            adata_s10_clean = remove_genes_by_plaque_correlation(
                adata_s10_qc, corr_s10, r2_threshold
            )
        else:
            adata_s10_clean = adata_s10_qc
    else:
        adata_s10_clean = adata_s10_qc

    if corr_s18 is not None:
        neg_probes_df = corr_s18[corr_s18['gene'].str.startswith('NegPrb')]
        if len(neg_probes_df) >= 2:
            r2_threshold = neg_probes_df.iloc[1]['R2']
            print(f"\nS18 R² threshold (from 2nd highest neg probe): {r2_threshold:.3f}")
            adata_s18_clean = remove_genes_by_plaque_correlation(
                adata_s18_qc, corr_s18, r2_threshold
            )
        else:
            adata_s18_clean = adata_s18_qc
    else:
        adata_s18_clean = adata_s18_qc

    print("\n" + "=" * 80)
    print("STEP 5: Combine sections")
    print("=" * 80)

    # Combine both sections
    adata_combined = ad.concat(
        [adata_s10_clean, adata_s18_clean],
        label='section',
        keys=['S10', 'S18'],
        index_unique=None
    )

    print(f"\nFinal combined data:")
    print(f"Cells: {adata_combined.n_obs}")
    print(f"Genes: {adata_combined.n_vars}")
    print(f"\nTarget from paper: 37,840 cells")
    print(f"Your result: {adata_combined.n_obs} cells")
    print(f"Difference: {37840 - adata_combined.n_obs} cells")

    return adata_combined, adata_s10_clean, adata_s18_clean


# Run the complete workflow
if __name__ == "__main__":
    DATA_PATH = '/path/to/your/data'  # Update this path

    adata_all, adata_s10, adata_s18 = complete_qc_workflow(DATA_PATH)

    # Note: The paper mentions additional microglial subclustering
    # to remove 141 contaminated cells. This would be done AFTER
    # the steps above using clustering analysis on microglial cells only.