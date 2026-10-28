"""Metacell-level CNV heatmaps (rows = metacells, columns = genomic bins)."""
import numpy as np, pandas as pd, matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, ListedColormap
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import linkage, leaves_list
from .plotting import _chr_ticks

EPI = ['luminal', 'basal', 'club', 'hillock']
STRO = ['endothelial', 'fibroblast', 'smooth_muscle', 'schwann']
SUB_COL = {'luminal': '#e08214', 'basal': '#8073ac', 'club': '#35978f', 'hillock': '#c51b7d',
           'stroma (reference)': '#4d4d4d', 'immune (reference)': '#bdbdbd'}

def chr_ticks(bins):
    return _chr_ticks(bins)

def _order(M):
    if len(M) < 3: return np.arange(len(M))
    return leaves_list(linkage(M, 'ward'))

def wong_heatmap(outs, bins, variant, which, fname, max_epi=45, max_ref=8, vmax=0.6, seed=0, title=None):
    rng = np.random.default_rng(seed); pats = sorted(outs)
    blocks = []
    for p in pats:
        d = outs[p]['variants'][variant]; mt = d['meta_table'].copy(); M = d[which]
        ct = mt.group.str.split('|').str[0].values
        # fraction of each metacell's cells coming from a tumour-labelled library (annotation only; never used in the analysis)
        lib = outs[p]['lib']; cm = d['meta']; fT = np.full(len(mt), np.nan)
        ok = cm >= 0; s = pd.Series((lib[ok] == 'T').astype(float)).groupby(cm[ok]).mean(); fT[s.index.values] = s.values
        epi = np.where(np.isin(ct, EPI) & ~mt.is_ref.values & ~mt.group.str.contains('_ho').values)[0]
        if len(epi) > max_epi: epi = np.sort(rng.choice(epi, max_epi, replace=False))
        epi = epi[_order(M[epi])]
        ref_s = np.where(mt.is_ref.values & np.isin(ct, STRO))[0]; ref_i = np.where(mt.is_ref.values & ~np.isin(ct, STRO))[0]
        ref_s = np.sort(rng.choice(ref_s, min(max_ref, len(ref_s)), replace=False)) if len(ref_s) else ref_s
        ref_i = np.sort(rng.choice(ref_i, min(max_ref, len(ref_i)), replace=False)) if len(ref_i) else ref_i
        rows = np.r_[epi, ref_s, ref_i].astype(int)
        lab = [ct[i] for i in epi] + ['stroma (reference)'] * len(ref_s) + ['immune (reference)'] * len(ref_i)
        blocks.append(dict(p=p, rows=rows, lab=lab, fT=fT[rows], M=M[rows], n_epi=len(epi), n_all_epi=int(np.isin(ct, EPI).sum() - mt.group.str.contains('_ho').values[np.isin(ct, EPI)].sum())))
    heights = [len(b['rows']) + 3 for b in blocks]
    H = 0.019 * sum(heights) + 1.9; fig = plt.figure(figsize=(7.2, H))
    gs = fig.add_gridspec(len(blocks), 3, height_ratios=heights, width_ratios=[0.012, 0.012, 1], hspace=0.12, wspace=0.01,
                          left=0.085, right=0.915, top=0.955, bottom=1.15 / H)
    bnd, mid, labx = chr_ticks(bins); norm = TwoSlopeNorm(0, -vmax, vmax)
    gray = mpl.colormaps['Greys']
    for k, b in enumerate(blocks):
        a0 = fig.add_subplot(gs[k, 0]); a1 = fig.add_subplot(gs[k, 1]); ax = fig.add_subplot(gs[k, 2])
        im = ax.imshow(b['M'], aspect='auto', cmap='RdBu_r', norm=norm, interpolation='nearest')
        a0.imshow(np.array([[mpl.colors.to_rgb(SUB_COL[l])] for l in b['lab']]), aspect='auto', interpolation='nearest')
        fT = np.where(np.isnan(b['fT']), 0.5, b['fT']); isref = np.array(['reference' in l for l in b['lab']])
        rgb = np.array([[gray(0.08 + 0.85 * f)[:3]] for f in fT]); rgb[isref] = 1.0
        a1.imshow(rgb, aspect='auto', interpolation='nearest')
        for a in (a0, a1): a.set_xticks([]); a.set_yticks([])
        for b_ in bnd: ax.axvline(b_, color='k', lw=0.25, alpha=0.5)
        n = b['n_epi']; ax.axhline(n - 0.5, color='k', lw=0.7)
        ax.set_yticks([]); ax.set_xticks(mid if k == len(blocks) - 1 else [])
        if k == len(blocks) - 1: ax.set_xticklabels(labx, fontsize=6); ax.tick_params(axis='x', length=0)
        ax.text(1.005, 1 - 0.5 * n / len(b['rows']), f"{b['p']}\n{b['n_epi']} of {b['n_all_epi']}\nepithelial\nmetacells", transform=ax.transAxes, fontsize=6, va='center', ha='left')
    cax = fig.add_axes([0.70, 0.62 / H, 0.215, 0.09 / H]); cb = fig.colorbar(im, cax=cax, orientation='horizontal'); cb.ax.tick_params(labelsize=6, length=2)
    cb.ax.set_title('log2 ratio vs stroma / immune reference', fontsize=6, pad=3)
    handles = [Patch(color=SUB_COL[k], label=k) for k in SUB_COL]
    handles += [Patch(facecolor=gray(0.93), label='cells from tumour-labelled library (annotation only)'), Patch(facecolor=gray(0.08), edgecolor='0.6', label='cells from benign-labelled library (annotation only)')]
    fig.legend(handles=handles, loc='lower left', ncol=3, fontsize=5.5, frameon=False, bbox_to_anchor=(0.075, 0.0), columnspacing=1.0, handlelength=1.0)
    if title: fig.suptitle(title, x=0.085, ha='left', fontsize=8, y=0.985)
    fig.savefig(fname, dpi=200); return fig, blocks
