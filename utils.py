from anndata import AnnData


def get_adata_section_id(adata: AnnData) -> str:
    return adata.uns.get('section_id', 'unknown')
