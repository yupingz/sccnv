"""Small simulator used by the unit tests and the quickstart example: a toy genome and Poisson counts in which 'tumour' cells carry
a gain and a loss on chromosome arms relative to diploid reference cells."""
import tempfile, pathlib
import numpy as np, pandas as pd, scipy.sparse as sp
from .core import load_gene_table

CHROMS = ['1', '2', '3', '4']


def toy_genome(genes_per_chr=800):
    """4 chromosomes x genes_per_chr genes, centromere at the midpoint. Returns (gene_table, centromeres, symbols, chromosome_of_gene)."""
    rows = []
    for c in CHROMS:
        for i in range(genes_per_chr):
            s = 1_000_000 + i * 100_000
            rows.append((f'G{c}_{i}', f'G{c}_{i}', c, s, s + 5000, '+', 'protein_coding'))
    t = pd.DataFrame(rows, columns=['gene_id', 'symbol', 'chr', 'start', 'end', 'strand', 'biotype'])
    p = pathlib.Path(tempfile.mkdtemp()) / 'genes.csv'; t.to_csv(p, index=False)
    mid = 1_000_000 + (genes_per_chr // 2) * 100_000
    return load_gene_table(p), {c: (mid, mid + 1) for c in CHROMS}, list(t.symbol), t.chr.values


def simulate(genome, n_ref=420, n_tum=240, depth=2500, gain=('2', 'q', 1.6), loss=('3', 'p', 0.5), seed=0):
    """Counts (cells x genes, CSR) with n_ref diploid reference cells followed by n_tum tumour cells.
    gain / loss = (chromosome, arm, fold change). 25% of the reference cells are flagged held-out (labelled '..._ho')."""
    table, cen, syms, chr_of = genome
    rng = np.random.default_rng(seed); n_g = len(syms)
    mu = rng.lognormal(0, 1.0, n_g); mu /= mu.sum()
    pos = table.loc[syms, 'pos'].values; arm = np.where(pos <= cen[chr_of[0]][0], 'p', 'q')
    f = np.ones(n_g)
    f[(chr_of == gain[0]) & (arm == gain[1])] = gain[2]
    f[(chr_of == loss[0]) & (arm == loss[1])] = loss[2]

    def cells(n, fac):
        p = mu * fac; p = p / p.sum()
        return rng.poisson(rng.gamma(20, depth / 20, n)[:, None] * p[None])

    X = sp.csr_matrix(np.vstack([cells(n_ref, np.ones(n_g)), cells(n_tum, f)]).astype(np.float32))
    kind = np.array(['ref'] * n_ref + ['tum'] * n_tum)
    ho = np.zeros(len(kind), bool); ho[rng.choice(n_ref, int(0.25 * n_ref), replace=False)] = True
    grp = np.array([f'{k}|c0' + ('_ho' if h else '') for k, h in zip(kind, ho)], dtype=object)
    return dict(X=X, genes=syms, kind=kind, ho=ho, grp=grp, ref_mask=(kind == 'ref') & ~ho, f=f, arm=arm, chr=chr_of)
