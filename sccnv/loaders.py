"""Load Seurat exports (counts CSC pieces written from R) into cells x genes CSR + metadata."""
import numpy as np, pandas as pd, scipy.sparse as sp, os

def load_export(d, x_file=None):
    d = d.rstrip('/')
    genes = [l.rstrip('\n') for l in open(f'{d}/genes.txt')]; cells = [l.rstrip('\n') for l in open(f'{d}/cells.txt')]
    i = np.fromfile(f'{d}/counts_i.bin', dtype=np.int32); p = np.fromfile(f'{d}/counts_p.bin', dtype=np.int32)
    if x_file is None:
        x_file = [f for f in os.listdir(d) if f.startswith('counts_x')][0]
    xf = f'{d}/{x_file}'
    x = np.fromfile(xf, dtype=np.int32 if 'int32' in x_file else np.float64).astype(np.float32)
    X = sp.csc_matrix((x, i, p), shape=(len(genes), len(cells))).T.tocsr()
    meta = pd.read_csv(f'{d}/meta.csv').set_index('cell').loc[cells]
    return X, genes, cells, meta
