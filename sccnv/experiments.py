"""Depth-downsampling experiment: does metacell aggregation preserve CNV sensitivity at low UMI depth?"""
import time, numpy as np, pandas as pd
from . import core as sc_cnv, benchmark as bench
from .core import Config, run, downsample_counts, infercnv_like

DEFAULT_METHODS = {
    'single-cell (pipeline)':        dict(kind='pipeline', cfg=dict(metacell_cells=1, thin=False)),
    'metacell-20 (pipeline)':        dict(kind='pipeline', cfg=dict(metacell_cells=20, thin=True)),
    'metacell-50 (pipeline)':        dict(kind='pipeline', cfg=dict(metacell_cells=50, thin=True)),
    'metacell-20, no count split':   dict(kind='pipeline', cfg=dict(metacell_cells=20, thin=False)),
    'inferCNV-style (single cell)':  dict(kind='baseline'),
}

def run_depth_experiment(X, genes, table, cen, grp, ref_masks, sel, cell_key, truth_by_key, norm_mask, fracs,
                         methods=None, base_cfg=None, seed=0, profiles=('seg', 'r'), log=print, run_kw=None, eval_kw=None):
    methods = methods or DEFAULT_METHODS; base_cfg = base_cfg or {}
    rng = np.random.default_rng(seed); rows = []
    ref_union = np.any([np.asarray(m) for m in ref_masks.values()], 0)
    for f in fracs:
        Xf = downsample_counts(X, f, rng)
        med = float(np.median(np.asarray(Xf.sum(1)).ravel()))
        for name, spec in methods.items():
            t0 = time.time()
            if spec['kind'] == 'baseline':
                P = {'baseline': infercnv_like(Xf, sel, ref_union)}
            else:
                cfg = Config(**{**base_cfg, **spec['cfg']})
                res = run(Xf, genes, table, cen, grp, ref_masks, cfg, verbose=False, sel=sel, **(run_kw or {}))
                P = {p: res.cell_matrix(p) for p in profiles}
            for pname, M in P.items():
                for d in bench.eval_profiles(M, cell_key, truth_by_key, norm_mask, **(eval_kw or {})):
                    rows.append(dict(frac=f, median_umi=med, method=name, profile=pname, **d))
            log(f"frac={f} umi~{med:.0f} {name}: {time.time()-t0:.0f}s")
    return pd.DataFrame(rows)
