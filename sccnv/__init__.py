"""sccnv: expression-only copy-number inference from scRNA-seq with lineage-aware references and count-split metacells."""
from .core import Config, run, select_genes_and_bins, load_gene_table, load_centromeres, cluster_cells
__all__ = ['Config', 'run', 'select_genes_and_bins', 'load_gene_table', 'load_centromeres', 'cluster_cells']
__version__ = '0.1.0'
