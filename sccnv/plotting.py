"""Plot helpers (matplotlib)."""
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

def _chr_ticks(bins, min_gap_frac=0.024):
    """chromosome boundaries, label positions and labels; labels of neighbouring chromosomes closer than min_gap_frac of the axis are merged ('20-22')."""
    chr_i = bins.chr_i.values
    bnd = np.where(np.diff(chr_i) != 0)[0] + 0.5
    mid = bins.groupby('chr_i').apply(lambda d: (d.index.min() + d.index.max()) / 2, include_groups=False).values
    lab = bins.groupby('chr_i').chr.first().astype(str).values
    gap = min_gap_frac * len(bins); groups = [[0]]
    for k in range(1, len(mid)):
        (groups[-1].append(k) if mid[k] - mid[groups[-1][-1]] < gap else groups.append([k]))
    pos = np.array([np.mean(mid[g]) for g in groups])
    txt = np.array([lab[g[0]] if len(g) == 1 else f"{lab[g[0]]}\u2013{lab[g[-1]]}" for g in groups], dtype=object)
    return bnd, pos, txt

def heatmap(res, which='seg', group_order=None, max_rows=30, vmax=0.6, ax=None, title='', seed=0, label_map=None):
    """metacells x bins heatmap, grouped by `group` (subsampled to max_rows per group)."""
    rng = np.random.default_rng(seed); mt = res.meta_table; M = getattr(res, which)
    groups = group_order or list(pd.unique(mt.group))
    rows, labs = [], []
    for g in groups:
        idx = np.where(mt.group.values == g)[0]
        if len(idx) == 0: continue
        if len(idx) > max_rows: idx = np.sort(rng.choice(idx, max_rows, replace=False))
        # order by similarity: first PC of profiles
        if len(idx) > 2:
            u = np.linalg.svd(M[idx] - M[idx].mean(0), full_matrices=False)[0][:, 0]; idx = idx[np.argsort(u)]
        rows.append(idx); labs.append((g, len(idx)))
    order = np.concatenate(rows); img = M[order]
    if ax is None: fig, ax = plt.subplots(figsize=(11, 6))
    ax.imshow(img, aspect='auto', cmap='RdBu_r', norm=TwoSlopeNorm(0, -vmax, vmax), interpolation='nearest')
    bnd, mid, lab = _chr_ticks(res.bins)
    for b in bnd: ax.axvline(b, color='k', lw=.4, alpha=.5)
    ax.set_xticks(mid); ax.set_xticklabels(lab, fontsize=7)
    pos = 0
    for g, n in labs:
        ax.axhline(pos - .5, color='k', lw=.6); ax.text(-3, pos + n / 2, (label_map or {}).get(g, g), ha='right', va='center', fontsize=7); pos += n
    ax.set_yticks([]); ax.set_title(title, fontsize=9)
    return ax
