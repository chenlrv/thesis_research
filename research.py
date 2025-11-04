import anndata as ad
import squidpy as sq
from anndata import AnnData

from qc import qc, remove_plaque_correlated_genes, remove_negative_probes, get_negative_probes

DATA_PATH_TEMPLATE = 'D:/thesis_research/resources/{section_id}'


def load_adata(section_id: str) -> AnnData:
    print(f'Loading adata for section {section_id}')

    adata = sq.read.nanostring(
        path=DATA_PATH_TEMPLATE.format(section_id=section_id),
        counts_file=f'exprMat_file.csv',
        meta_file=f'metadata_file.csv',
        fov_file=f'fov_positions_file.csv')
    adata.uns['section_id'] = section_id

    print(f'Adata {section_id} initial shapes:\n'
          f'{adata}\n'
          f'{adata.var.head()}\n'
          f'Cells: {adata.n_obs}\n'
          f'Genes: {adata.n_vars}')
    return adata


def process_adata(section_id: str) -> tuple[AnnData, AnnData]:
    adata = load_adata(section_id)
    negative_probes = get_negative_probes(adata)
    adata = remove_negative_probes(adata, negative_probes)
    adata = qc(adata)
    return adata, negative_probes


def run():
    s10_adata, s10_negative_probes = process_adata(section_id='GSM8199188_ID61-ID62_S10')
    s18_adata, s18_negative_probes = process_adata(section_id='GSM8199189_ID67-ID68_S18')

    sections_adata_unified = s10_adata.concatenate(
        s18_adata,
        batch_key='section_id',
        batch_categories=['S10', 'S18'],
        index_unique='-')

    negative_probes_unified = s10_negative_probes.concatenate(
        s18_negative_probes,
        batch_key='section_id',
        batch_categories=['S10', 'S18'],
        index_unique='-')

    final_adata = ad.concat(
        [sections_adata_unified, negative_probes_unified],
        axis=1,
        merge='same'
    )
    filtered_adata = remove_plaque_correlated_genes(final_adata)


run()
