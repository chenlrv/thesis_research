import pandas as pd
import json


CELL_TYPE_NAMES = {
    'VLMC': 'VLMC',
    'Peri': 'Pericyte',
    'SMC': 'SMC',
    'Endo': 'Endothelial',
    'Microglia': 'Microglia',
    'BAM': 'BAM',
    'Astro': 'Astrocyte',
    'Oligo': 'Oligodendrocyte',
    'Glut': 'Glutamatergic',
    'Gaba': 'GABAergic',
    'OPC': 'OPC',
    'Epen': 'Ependymal',
    'Macro': 'Macrophage',
    'Mono': 'Monocyte',
    'DC': 'Dendritic Cell',
}


def get_allen_brain_atlas_cell_type_markers(excel_file_path: str = None):
    """
    Download and process Allen Brain Atlas cell type markers
    """

    if excel_file_path:
        print(f'Loading data from {excel_file_path}...')
        df = pd.read_excel(excel_file_path)
    else:
        url = "https://allen-brain-cell-atlas.s3-us-west-2.amazonaws.com/metadata/WMB-taxonomy/20230830/cl.df_CCN202307220.xlsx"
        df = pd.read_excel(url)

        df.to_csv('allen_markers.xlsx', index=False)

    print(f'Loaded {len(df)} cell clusters from Allen Brain Atlas.')

    unknown_cell_types = get_unknown_cell_types(df)
    add_cell_type(df, unknown_cell_types)
    cell_type_markers = extract_cell_type_markers(df)
    print('Extracted markers for cell types:', list(cell_type_markers.keys()))


def get_unknown_cell_types(df: pd.DataFrame) -> set[str]:
    unknown_types = set()

    for idx, row in df.iterrows():
        supertype_label = row.get('supertype_label', None)
        cell_type = extract_cell_type(supertype_label)
        if cell_type == supertype_label:
            unknown_types.add(supertype_label)

    if unknown_types:
        print('Warning: Unrecognized supertype_labels found:', unknown_types)

    return unknown_types

def add_cell_type(df: pd.DataFrame, unknown_cell_types: set[str]) -> None:
    for idx, row in df.iterrows():
        supertype_label = row.get('supertype_label', None)
        if supertype_label in unknown_cell_types:
            df.at[idx, 'cell_type'] = None
        else:
            df.at[idx, 'cell_type'] = extract_cell_type(supertype_label)


def extract_cell_type(supertype_label: str) -> str:
    """
    Extract cell type from supertype_label

    'VLMC NN_2' -> 'VLMC'
    'Microglia NN_1' -> 'Microglia'
    """

    for key in CELL_TYPE_NAMES.keys():
        if key.lower() in str(supertype_label).lower():
            return CELL_TYPE_NAMES[key]

    return supertype_label


def extract_cell_type_markers(df: pd.DataFrame) -> dict[str, list]:
    cell_type_markers = {}

    for cell_type in df['cell_type'].dropna().unique():
        type_df = df[df['cell_type'] == cell_type]

        all_markers = set()
        for markers_str in type_df['cluster.markers.combo']:
            if pd.notna(markers_str):
                all_markers.update(markers_str.split())

        cell_type_markers[cell_type] = sorted(list(all_markers))

    return cell_type_markers


# def save_outputs(df):
#     """
#     Save organized outputs
#     """
#
#     print("\nCreating outputs...")
#
#     # 1. Simple CSV: cell_type | cell_type_id | label | markers
#     output = df[[
#         'cell_type',
#         'cl',
#         'cluster_id_label',
#         'cluster.markers.combo'
#     ]].copy()
#
#     output.columns = ['cell_type', 'cell_type_id', 'label', 'markers']
#     output = output.sort_values('cell_type')
#     output.to_csv('cell_types_simple.csv', index=False)
#     print("✓ Saved: cell_types_simple.csv")
#
#     # 2. Dictionary: cell_type -> markers
#     markers_dict = {}
#     for cell_type in df['cell_type'].unique():
#         type_df = df[df['cell_type'] == cell_type]
#
#         # Get all unique markers
#         all_markers = set()
#         for markers_str in type_df['cluster.markers.combo']:
#             if pd.notna(markers_str):
#                 all_markers.update(markers_str.split())
#
#         markers_dict[cell_type] = sorted(list(all_markers))
#
#     with open('cell_type_markers.json', 'w') as f:
#         json.dump(markers_dict, f, indent=2)
#     print("✓ Saved: cell_type_markers.json")
#
#     return markers_dict
#
#
# def main(excel_file_path):
#     """
#     Main function
#     """
#
#     df = categorize_cells(excel_file_path)
#     markers = save_outputs(df)
#
#     # Show examples
#     print("\n" + "=" * 70)
#     print("EXAMPLES")
#     print("=" * 70)
#
#     for cell_type in ['Oligodendrocyte', 'Microglia', 'Astrocyte', 'Pericyte', 'Endothelial']:
#         if cell_type in df['cell_type'].values:
#             type_df = df[df['cell_type'] == cell_type]
#             print(f"\n{cell_type}: {len(type_df)} subtypes")
#             print(f"  Example IDs: {type_df['cl'].head(3).tolist()}")
#             if cell_type in markers:
#                 print(f"  Markers ({len(markers[cell_type])}): {markers[cell_type][:10]}")
#
#     print("\n" + "=" * 70)
#     print("USAGE")
#     print("=" * 70)




get_allen_brain_atlas_cell_type_markers()